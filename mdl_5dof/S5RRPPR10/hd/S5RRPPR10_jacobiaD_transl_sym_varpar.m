% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPR10
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPPR10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPR10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:31
	% EndTime: 2019-12-31 19:45:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:32
	% EndTime: 2019-12-31 19:45:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:32
	% EndTime: 2019-12-31 19:45:32
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(6) + r_i_i_C(3);
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
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:32
	% EndTime: 2019-12-31 19:45:32
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (55->22), mult. (178->38), div. (0->0), fcn. (133->6), ass. (0->17)
	t160 = sin(qJ(2));
	t162 = cos(qJ(2));
	t158 = sin(pkin(8));
	t159 = cos(pkin(8));
	t169 = r_i_i_C(1) * t159 - r_i_i_C(2) * t158 + pkin(2);
	t174 = r_i_i_C(3) + qJ(3);
	t175 = t169 * t160 - t174 * t162;
	t176 = t175 * qJD(2) - t160 * qJD(3);
	t161 = sin(qJ(1));
	t173 = qJD(1) * t161;
	t163 = cos(qJ(1));
	t172 = qJD(1) * t163;
	t171 = qJD(2) * t163;
	t168 = t158 * r_i_i_C(1) + t159 * r_i_i_C(2) + pkin(6);
	t167 = -t174 * t160 - t169 * t162;
	t165 = -pkin(1) + t167;
	t1 = [t176 * t161 + (-t168 * t161 + t165 * t163) * qJD(1), (t169 * t173 - t174 * t171) * t160 + (-t174 * t173 + (-t169 * qJD(2) + qJD(3)) * t163) * t162, -t160 * t173 + t162 * t171, 0, 0; -t176 * t163 + (t165 * t161 + t168 * t163) * qJD(1), -t175 * t172 + (t167 * qJD(2) + qJD(3) * t162) * t161, t161 * qJD(2) * t162 + t160 * t172, 0, 0; 0, -t176, qJD(2) * t160, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:32
	% EndTime: 2019-12-31 19:45:33
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (88->38), mult. (286->64), div. (0->0), fcn. (229->6), ass. (0->32)
	t179 = sin(qJ(2));
	t177 = sin(pkin(8));
	t178 = cos(pkin(8));
	t203 = r_i_i_C(3) + qJ(4);
	t206 = pkin(3) + r_i_i_C(1);
	t207 = t203 * t177 + t206 * t178 + pkin(2);
	t181 = cos(qJ(2));
	t204 = r_i_i_C(2) + qJ(3);
	t208 = t204 * t181;
	t211 = t207 * t179 - t208;
	t194 = qJD(4) * t177;
	t187 = t179 * qJD(3) + t181 * t194;
	t210 = (-pkin(2) * t179 + t208) * qJD(2) + t187;
	t182 = cos(qJ(1));
	t202 = t177 * t182;
	t180 = sin(qJ(1));
	t201 = t180 * t181;
	t200 = t182 * t178;
	t199 = qJD(1) * t180;
	t198 = qJD(1) * t182;
	t197 = qJD(2) * t179;
	t196 = qJD(2) * t181;
	t195 = qJD(2) * t182;
	t193 = t178 * qJD(4);
	t192 = t180 * t197;
	t191 = t179 * t195;
	t190 = t204 * t179;
	t186 = -pkin(2) * t181 - pkin(1) - t190;
	t183 = -t179 * t194 + qJD(3) * t181 + (-t181 * t207 - t190) * qJD(2);
	t175 = -t177 * t192 + (-t180 * t178 + t181 * t202) * qJD(1);
	t173 = t177 * t191 + (t177 * t201 + t200) * qJD(1);
	t1 = [-t182 * t193 + t206 * (t178 * t192 + (-t177 * t180 - t181 * t200) * qJD(1)) - t203 * t175 - t210 * t180 + (-t180 * pkin(6) + t186 * t182) * qJD(1), t183 * t182 + t211 * t199, -t179 * t199 + t181 * t195, -t173, 0; -t180 * t193 + t206 * (-t178 * t191 + (-t178 * t201 + t202) * qJD(1)) - t203 * t173 + t210 * t182 + (t182 * pkin(6) + t186 * t180) * qJD(1), t183 * t180 - t198 * t211, t179 * t198 + t180 * t196, t175, 0; 0, -qJD(2) * t211 + t187, t197, t177 * t196, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:33
	% EndTime: 2019-12-31 19:45:33
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (203->66), mult. (659->110), div. (0->0), fcn. (600->8), ass. (0->48)
	t241 = sin(qJ(2));
	t238 = sin(pkin(8));
	t239 = cos(pkin(8));
	t240 = sin(qJ(5));
	t243 = cos(qJ(5));
	t275 = pkin(3) + pkin(4);
	t251 = t243 * r_i_i_C(1) - t240 * r_i_i_C(2) + t275;
	t252 = t240 * r_i_i_C(1) + t243 * r_i_i_C(2) + qJ(4);
	t276 = t252 * t238 + t251 * t239 + pkin(2);
	t244 = cos(qJ(2));
	t263 = -r_i_i_C(3) - pkin(7) + qJ(3);
	t278 = t263 * t244;
	t281 = t276 * t241 - t278;
	t264 = t241 * qJD(3);
	t280 = (-t241 * pkin(2) + t278) * qJD(2) + t264;
	t253 = t238 * t240 + t239 * t243;
	t254 = t238 * t243 - t239 * t240;
	t249 = t254 * r_i_i_C(1) - t253 * r_i_i_C(2);
	t277 = t238 * qJD(4) + t249 * qJD(5);
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
	t262 = t244 * t272;
	t261 = t242 * t268;
	t260 = t241 * t266;
	t259 = t263 * t241;
	t230 = t238 * t273 + t271;
	t231 = t239 * t273 - t272;
	t256 = -t230 * t243 + t231 * t240;
	t255 = t230 * t240 + t231 * t243;
	t233 = t242 * t238 + t244 * t271;
	t250 = -pkin(2) * t244 - pkin(1) - t259;
	t246 = t244 * qJD(3) - t277 * t241 + (-t244 * t276 - t259) * qJD(2);
	t232 = -t242 * t239 + t262;
	t229 = t233 * qJD(1) - t239 * t261;
	t228 = qJD(1) * t262 - t238 * t261 - t239 * t270;
	t227 = -t231 * qJD(1) - t239 * t260;
	t226 = t230 * qJD(1) + t238 * t260;
	t225 = -t226 * t240 + t227 * t243 + (t232 * t243 - t233 * t240) * qJD(5);
	t224 = -t226 * t243 - t227 * t240 + (-t232 * t240 - t233 * t243) * qJD(5);
	t1 = [-t230 * qJD(4) - t252 * t228 - t251 * t229 + (t256 * r_i_i_C(1) + t255 * r_i_i_C(2)) * qJD(5) - t280 * t242 + (-t242 * pkin(6) + t250 * t245) * qJD(1), t246 * t245 + t281 * t270, -t241 * t270 + t244 * t266, -t226, t224 * r_i_i_C(1) - t225 * r_i_i_C(2); t225 * r_i_i_C(1) + t224 * r_i_i_C(2) - t226 * qJ(4) + t232 * qJD(4) + t275 * t227 + t280 * t245 + (pkin(6) * t245 + t250 * t242) * qJD(1), t246 * t242 - t269 * t281, t241 * t269 + t242 * t267, t228, (t228 * t243 - t229 * t240) * r_i_i_C(1) + (-t228 * t240 - t229 * t243) * r_i_i_C(2) + (-t255 * r_i_i_C(1) + t256 * r_i_i_C(2)) * qJD(5); 0, -qJD(2) * t281 + t277 * t244 + t264, t268, t238 * t267, (-t253 * r_i_i_C(1) - t254 * r_i_i_C(2)) * t241 * qJD(5) + t249 * t267;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end