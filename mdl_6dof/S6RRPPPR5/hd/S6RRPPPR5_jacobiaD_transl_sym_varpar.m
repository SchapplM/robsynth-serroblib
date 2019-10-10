% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPPR5
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:25
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPPR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(7) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0, 0; 0, -t21, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:10
	% EndTime: 2019-10-10 09:25:10
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (55->22), mult. (178->38), div. (0->0), fcn. (133->6), ass. (0->17)
	t160 = sin(qJ(2));
	t162 = cos(qJ(2));
	t158 = sin(pkin(9));
	t159 = cos(pkin(9));
	t169 = r_i_i_C(1) * t159 - r_i_i_C(2) * t158 + pkin(2);
	t174 = r_i_i_C(3) + qJ(3);
	t175 = t169 * t160 - t174 * t162;
	t176 = t175 * qJD(2) - t160 * qJD(3);
	t161 = sin(qJ(1));
	t173 = qJD(1) * t161;
	t163 = cos(qJ(1));
	t172 = qJD(1) * t163;
	t171 = qJD(2) * t163;
	t168 = t158 * r_i_i_C(1) + t159 * r_i_i_C(2) + pkin(7);
	t167 = -t174 * t160 - t169 * t162;
	t165 = -pkin(1) + t167;
	t1 = [t176 * t161 + (-t168 * t161 + t165 * t163) * qJD(1), (t169 * t173 - t174 * t171) * t160 + (-t174 * t173 + (-t169 * qJD(2) + qJD(3)) * t163) * t162, -t160 * t173 + t162 * t171, 0, 0, 0; -t176 * t163 + (t165 * t161 + t168 * t163) * qJD(1), -t175 * t172 + (t167 * qJD(2) + qJD(3) * t162) * t161, t161 * qJD(2) * t162 + t160 * t172, 0, 0, 0; 0, -t176, qJD(2) * t160, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:10
	% EndTime: 2019-10-10 09:25:10
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (88->38), mult. (286->64), div. (0->0), fcn. (229->6), ass. (0->32)
	t184 = sin(qJ(2));
	t182 = sin(pkin(9));
	t183 = cos(pkin(9));
	t208 = r_i_i_C(3) + qJ(4);
	t211 = pkin(3) - r_i_i_C(2);
	t212 = t208 * t182 + t211 * t183 + pkin(2);
	t186 = cos(qJ(2));
	t209 = r_i_i_C(1) + qJ(3);
	t213 = t209 * t186;
	t216 = t212 * t184 - t213;
	t199 = qJD(4) * t182;
	t192 = t184 * qJD(3) + t186 * t199;
	t215 = (-pkin(2) * t184 + t213) * qJD(2) + t192;
	t187 = cos(qJ(1));
	t207 = t182 * t187;
	t185 = sin(qJ(1));
	t206 = t185 * t186;
	t205 = t187 * t183;
	t204 = qJD(1) * t185;
	t203 = qJD(1) * t187;
	t202 = qJD(2) * t184;
	t201 = qJD(2) * t186;
	t200 = qJD(2) * t187;
	t198 = t183 * qJD(4);
	t197 = t185 * t202;
	t196 = t184 * t200;
	t195 = t209 * t184;
	t191 = -pkin(2) * t186 - pkin(1) - t195;
	t188 = -t184 * t199 + qJD(3) * t186 + (-t186 * t212 - t195) * qJD(2);
	t180 = -t182 * t197 + (-t185 * t183 + t186 * t207) * qJD(1);
	t178 = t182 * t196 + (t182 * t206 + t205) * qJD(1);
	t1 = [-t187 * t198 - t211 * (-t183 * t197 + (t182 * t185 + t186 * t205) * qJD(1)) - t208 * t180 - t215 * t185 + (-t185 * pkin(7) + t191 * t187) * qJD(1), t188 * t187 + t216 * t204, -t184 * t204 + t186 * t200, -t178, 0, 0; -t185 * t198 - t211 * (t183 * t196 + (t183 * t206 - t207) * qJD(1)) - t208 * t178 + t215 * t187 + (t187 * pkin(7) + t191 * t185) * qJD(1), t188 * t185 - t203 * t216, t184 * t203 + t185 * t201, t180, 0, 0; 0, -qJD(2) * t216 + t192, t202, t182 * t201, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:10
	% EndTime: 2019-10-10 09:25:10
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (119->42), mult. (383->68), div. (0->0), fcn. (312->6), ass. (0->31)
	t180 = sin(qJ(2));
	t182 = cos(qJ(2));
	t194 = pkin(4) - r_i_i_C(2) + qJ(3);
	t178 = sin(pkin(9));
	t179 = cos(pkin(9));
	t195 = pkin(3) + r_i_i_C(3) + qJ(5);
	t203 = r_i_i_C(1) + qJ(4);
	t204 = t203 * t178 + t195 * t179 + pkin(2);
	t209 = t204 * t180 - t194 * t182;
	t188 = qJD(4) * t178 + qJD(5) * t179;
	t208 = (t194 * qJD(2) + t188) * t182 - (pkin(2) * qJD(2) - qJD(3)) * t180;
	t181 = sin(qJ(1));
	t202 = t181 * t178;
	t183 = cos(qJ(1));
	t201 = t182 * t183;
	t200 = qJD(1) * t181;
	t199 = qJD(1) * t183;
	t198 = qJD(2) * t180;
	t197 = qJD(2) * t182;
	t196 = qJD(2) * t183;
	t193 = t180 * t196;
	t192 = t181 * t198;
	t190 = t194 * t180;
	t189 = -t179 * qJD(4) + t178 * qJD(5);
	t187 = -pkin(2) * t182 - pkin(1) - t190;
	t184 = qJD(3) * t182 - t188 * t180 + (-t182 * t204 - t190) * qJD(2);
	t175 = -t179 * t192 + (t179 * t201 + t202) * qJD(1);
	t174 = -t178 * t192 + (t178 * t201 - t181 * t179) * qJD(1);
	t173 = -t178 * t199 + (t182 * t200 + t193) * t179;
	t172 = t178 * t193 + (t179 * t183 + t182 * t202) * qJD(1);
	t1 = [t189 * t183 - t203 * t174 - t195 * t175 - t208 * t181 + (-t181 * pkin(7) + t187 * t183) * qJD(1), t184 * t183 + t209 * t200, -t180 * t200 + t182 * t196, -t172, -t173, 0; t189 * t181 - t203 * t172 - t195 * t173 + t208 * t183 + (t183 * pkin(7) + t187 * t181) * qJD(1), t184 * t181 - t199 * t209, t180 * t199 + t181 * t197, t174, t175, 0; 0, -qJD(2) * t209 + t180 * qJD(3) + t188 * t182, t198, t178 * t197, t179 * t197, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:10
	% EndTime: 2019-10-10 09:25:11
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (234->72), mult. (756->115), div. (0->0), fcn. (683->8), ass. (0->49)
	t241 = sin(qJ(2));
	t238 = sin(pkin(9));
	t239 = cos(pkin(9));
	t240 = sin(qJ(6));
	t243 = cos(qJ(6));
	t274 = -pkin(5) - qJ(4);
	t251 = t243 * r_i_i_C(1) - t240 * r_i_i_C(2) - t274;
	t275 = -pkin(3) - qJ(5);
	t252 = t240 * r_i_i_C(1) + t243 * r_i_i_C(2) - t275;
	t277 = t251 * t238 + t252 * t239 + pkin(2);
	t244 = cos(qJ(2));
	t261 = r_i_i_C(3) + pkin(8) + pkin(4) + qJ(3);
	t279 = t261 * t244;
	t283 = t277 * t241 - t279;
	t265 = t241 * qJD(3);
	t282 = (-t241 * pkin(2) + t279) * qJD(2) + t265;
	t254 = t238 * t240 - t239 * t243;
	t255 = t238 * t243 + t239 * t240;
	t278 = t254 * r_i_i_C(1) + t255 * r_i_i_C(2);
	t281 = -t238 * qJD(4) - t239 * qJD(5) + t278 * qJD(6);
	t242 = sin(qJ(1));
	t273 = t242 * t244;
	t245 = cos(qJ(1));
	t272 = t245 * t238;
	t271 = t245 * t239;
	t270 = qJD(1) * t242;
	t269 = qJD(1) * t245;
	t268 = qJD(2) * t241;
	t267 = qJD(2) * t244;
	t266 = qJD(2) * t245;
	t264 = t238 * t269;
	t263 = t242 * t268;
	t262 = t241 * t266;
	t260 = t261 * t241;
	t229 = t238 * t273 + t271;
	t230 = t239 * t273 - t272;
	t257 = t229 * t243 + t230 * t240;
	t256 = t229 * t240 - t230 * t243;
	t232 = t242 * t238 + t244 * t271;
	t250 = -pkin(2) * t244 - pkin(1) - t260;
	t246 = t244 * qJD(3) + t281 * t241 + (-t244 * t277 - t260) * qJD(2);
	t231 = -t242 * t239 + t244 * t272;
	t228 = t232 * qJD(1) - t239 * t263;
	t227 = -t238 * t263 - t239 * t270 + t244 * t264;
	t226 = -t264 + (t244 * t270 + t262) * t239;
	t225 = t229 * qJD(1) + t238 * t262;
	t224 = -t225 * t243 - t226 * t240 + (-t231 * t240 + t232 * t243) * qJD(6);
	t223 = t225 * t240 - t226 * t243 + (-t231 * t243 - t232 * t240) * qJD(6);
	t1 = [-t229 * qJD(4) - t230 * qJD(5) - t252 * t228 - t251 * t227 + (t256 * r_i_i_C(1) + t257 * r_i_i_C(2)) * qJD(6) - t282 * t242 + (-t242 * pkin(7) + t250 * t245) * qJD(1), t246 * t245 + t283 * t270, -t241 * t270 + t244 * t266, -t225, -t226, t223 * r_i_i_C(1) - t224 * r_i_i_C(2); t224 * r_i_i_C(1) + t223 * r_i_i_C(2) + t231 * qJD(4) + t232 * qJD(5) + t275 * t226 + t274 * t225 + t282 * t245 + (pkin(7) * t245 + t250 * t242) * qJD(1), t246 * t242 - t269 * t283, t241 * t269 + t242 * t267, t227, t228, (-t227 * t240 + t228 * t243) * r_i_i_C(1) + (-t227 * t243 - t228 * t240) * r_i_i_C(2) + (-t257 * r_i_i_C(1) + t256 * r_i_i_C(2)) * qJD(6); 0, -qJD(2) * t283 - t281 * t244 + t265, t268, t238 * t267, t239 * t267, (-t255 * r_i_i_C(1) + t254 * r_i_i_C(2)) * t241 * qJD(6) - t278 * t267;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end