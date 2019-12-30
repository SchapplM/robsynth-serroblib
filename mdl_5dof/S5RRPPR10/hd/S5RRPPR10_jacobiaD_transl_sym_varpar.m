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
% Datum: 2019-12-29 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
	% StartTime: 2019-12-29 18:31:20
	% EndTime: 2019-12-29 18:31:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:31:20
	% EndTime: 2019-12-29 18:31:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:31:20
	% EndTime: 2019-12-29 18:31:20
	% DurationCPUTime: 0.15s
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
	t1 = [t23 * t25 + (-t18 * t28 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:31:27
	% EndTime: 2019-12-29 18:31:27
	% DurationCPUTime: 0.31s
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
	% StartTime: 2019-12-29 18:31:27
	% EndTime: 2019-12-29 18:31:27
	% DurationCPUTime: 0.34s
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
	% StartTime: 2019-12-29 18:31:23
	% EndTime: 2019-12-29 18:31:23
	% DurationCPUTime: 0.83s
	% Computational Cost: add. (203->66), mult. (659->110), div. (0->0), fcn. (600->8), ass. (0->48)
	t243 = sin(qJ(2));
	t240 = sin(pkin(8));
	t241 = cos(pkin(8));
	t242 = sin(qJ(5));
	t245 = cos(qJ(5));
	t277 = pkin(3) + pkin(4);
	t253 = t245 * r_i_i_C(1) - t242 * r_i_i_C(2) + t277;
	t254 = t242 * r_i_i_C(1) + t245 * r_i_i_C(2) + qJ(4);
	t278 = t254 * t240 + t253 * t241 + pkin(2);
	t246 = cos(qJ(2));
	t265 = -r_i_i_C(3) - pkin(7) + qJ(3);
	t280 = t265 * t246;
	t283 = t278 * t243 - t280;
	t266 = t243 * qJD(3);
	t282 = (-t243 * pkin(2) + t280) * qJD(2) + t266;
	t255 = t240 * t242 + t241 * t245;
	t256 = t240 * t245 - t241 * t242;
	t251 = t256 * r_i_i_C(1) - t255 * r_i_i_C(2);
	t279 = t240 * qJD(4) + t251 * qJD(5);
	t244 = sin(qJ(1));
	t275 = t244 * t246;
	t247 = cos(qJ(1));
	t274 = t247 * t240;
	t273 = t247 * t241;
	t272 = qJD(1) * t244;
	t271 = qJD(1) * t247;
	t270 = qJD(2) * t243;
	t269 = qJD(2) * t246;
	t268 = qJD(2) * t247;
	t264 = t246 * t274;
	t263 = t244 * t270;
	t262 = t243 * t268;
	t261 = t265 * t243;
	t232 = t240 * t275 + t273;
	t233 = t241 * t275 - t274;
	t258 = -t232 * t245 + t233 * t242;
	t257 = t232 * t242 + t233 * t245;
	t235 = t244 * t240 + t246 * t273;
	t252 = -pkin(2) * t246 - pkin(1) - t261;
	t248 = t246 * qJD(3) - t279 * t243 + (-t246 * t278 - t261) * qJD(2);
	t234 = -t244 * t241 + t264;
	t231 = t235 * qJD(1) - t241 * t263;
	t230 = qJD(1) * t264 - t240 * t263 - t241 * t272;
	t229 = -t233 * qJD(1) - t241 * t262;
	t228 = t232 * qJD(1) + t240 * t262;
	t227 = -t228 * t242 + t229 * t245 + (t234 * t245 - t235 * t242) * qJD(5);
	t226 = -t228 * t245 - t229 * t242 + (-t234 * t242 - t235 * t245) * qJD(5);
	t1 = [-t232 * qJD(4) - t254 * t230 - t253 * t231 + (t258 * r_i_i_C(1) + t257 * r_i_i_C(2)) * qJD(5) - t282 * t244 + (-t244 * pkin(6) + t252 * t247) * qJD(1), t248 * t247 + t283 * t272, -t243 * t272 + t246 * t268, -t228, r_i_i_C(1) * t226 - r_i_i_C(2) * t227; t227 * r_i_i_C(1) + t226 * r_i_i_C(2) - t228 * qJ(4) + t234 * qJD(4) + t277 * t229 + t282 * t247 + (pkin(6) * t247 + t252 * t244) * qJD(1), t248 * t244 - t271 * t283, t243 * t271 + t244 * t269, t230, (t230 * t245 - t231 * t242) * r_i_i_C(1) + (-t230 * t242 - t231 * t245) * r_i_i_C(2) + (-t257 * r_i_i_C(1) + t258 * r_i_i_C(2)) * qJD(5); 0, -qJD(2) * t283 + t279 * t246 + t266, t270, t240 * t269, (-t255 * r_i_i_C(1) - t256 * r_i_i_C(2)) * t243 * qJD(5) + t251 * t269;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end