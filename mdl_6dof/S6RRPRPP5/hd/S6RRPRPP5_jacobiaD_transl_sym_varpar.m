% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:02
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPP5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:02:41
	% EndTime: 2019-10-10 10:02:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:02:41
	% EndTime: 2019-10-10 10:02:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:02:41
	% EndTime: 2019-10-10 10:02:41
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
	% StartTime: 2019-10-10 10:02:42
	% EndTime: 2019-10-10 10:02:42
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (44->20), mult. (134->34), div. (0->0), fcn. (92->4), ass. (0->15)
	t139 = sin(qJ(2));
	t141 = cos(qJ(2));
	t151 = r_i_i_C(3) + qJ(3);
	t153 = pkin(2) - r_i_i_C(2);
	t154 = t153 * t139 - t151 * t141;
	t155 = t154 * qJD(2) - t139 * qJD(3);
	t152 = pkin(7) + r_i_i_C(1);
	t140 = sin(qJ(1));
	t150 = qJD(1) * t140;
	t142 = cos(qJ(1));
	t149 = qJD(1) * t142;
	t148 = qJD(2) * t142;
	t146 = -t151 * t139 - t153 * t141;
	t144 = -pkin(1) + t146;
	t1 = [t155 * t140 + (-t152 * t140 + t144 * t142) * qJD(1), (-t151 * t148 + t153 * t150) * t139 + (-t151 * t150 + (-t153 * qJD(2) + qJD(3)) * t142) * t141, -t139 * t150 + t141 * t148, 0, 0, 0; -t155 * t142 + (t144 * t140 + t152 * t142) * qJD(1), -t154 * t149 + (t146 * qJD(2) + qJD(3) * t141) * t140, t140 * qJD(2) * t141 + t139 * t149, 0, 0, 0; 0, -t155, qJD(2) * t139, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:02:42
	% EndTime: 2019-10-10 10:02:43
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (101->43), mult. (318->73), div. (0->0), fcn. (248->6), ass. (0->34)
	t199 = sin(qJ(2));
	t202 = cos(qJ(2));
	t216 = pkin(2) + pkin(8) + r_i_i_C(3);
	t213 = t216 * t199;
	t229 = (-qJ(3) * t202 + t213) * qJD(2) - t199 * qJD(3);
	t201 = cos(qJ(4));
	t211 = qJD(4) * t199 + qJD(1);
	t227 = t201 * t211;
	t198 = sin(qJ(4));
	t226 = t211 * t198;
	t225 = pkin(3) + pkin(7);
	t200 = sin(qJ(1));
	t223 = qJD(1) * t200;
	t203 = cos(qJ(1));
	t222 = qJD(1) * t203;
	t221 = qJD(2) * t199;
	t220 = qJD(2) * t202;
	t219 = qJD(2) * t203;
	t218 = qJD(4) * t202;
	t215 = t200 * t220;
	t214 = t202 * t219;
	t212 = t216 * t202;
	t210 = -qJD(1) * t199 - qJD(4);
	t209 = t210 * t203;
	t208 = r_i_i_C(1) * t198 + r_i_i_C(2) * t201 + qJ(3);
	t207 = -qJ(3) * t199 - pkin(1) - t212;
	t206 = qJD(3) + (r_i_i_C(1) * t201 - r_i_i_C(2) * t198) * qJD(4);
	t205 = t210 * t200 + t214;
	t204 = t208 * t202 - t213;
	t197 = t205 * t198 + t203 * t227;
	t196 = t205 * t201 - t203 * t226;
	t195 = -t200 * t227 + (t209 - t215) * t198;
	t194 = t201 * t209 + (-t201 * t220 + t226) * t200;
	t1 = [t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t229 * t200 + (-t225 * t200 + t207 * t203) * qJD(1), (-t208 * t219 + t216 * t223) * t199 + (-t208 * t223 + (-t216 * qJD(2) + t206) * t203) * t202, -t199 * t223 + t214, t196 * r_i_i_C(1) - t197 * r_i_i_C(2), 0, 0; t197 * r_i_i_C(1) + t196 * r_i_i_C(2) - t229 * t203 + (t207 * t200 + t225 * t203) * qJD(1), t204 * t222 + (t206 * t202 + (-t208 * t199 - t212) * qJD(2)) * t200, t199 * t222 + t215, -t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0, 0; 0, t204 * qJD(2) + t206 * t199, t221, (-t198 * t221 + t201 * t218) * r_i_i_C(2) + (t198 * t218 + t201 * t221) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:02:43
	% EndTime: 2019-10-10 10:02:43
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (182->54), mult. (558->87), div. (0->0), fcn. (464->6), ass. (0->39)
	t241 = sin(qJ(2));
	t244 = cos(qJ(2));
	t264 = qJD(2) * t244;
	t243 = cos(qJ(4));
	t261 = pkin(2) + pkin(8) + r_i_i_C(2);
	t276 = t261 * qJD(2) + qJD(5) * t243 - qJD(3);
	t278 = -qJ(3) * t264 + t276 * t241;
	t242 = sin(qJ(1));
	t258 = t242 * t264;
	t245 = cos(qJ(1));
	t266 = qJD(1) * t245;
	t277 = t241 * t266 + t258;
	t240 = sin(qJ(4));
	t272 = r_i_i_C(3) + qJ(5);
	t273 = pkin(4) + r_i_i_C(1);
	t275 = t273 * t240 - t272 * t243 + qJ(3);
	t274 = pkin(3) + pkin(7);
	t271 = t240 * t245;
	t270 = t242 * t240;
	t269 = t242 * t243;
	t268 = t243 * t245;
	t267 = qJD(1) * t242;
	t265 = qJD(2) * t241;
	t263 = qJD(4) * t244;
	t262 = qJD(5) * t240;
	t257 = t245 * t264;
	t256 = qJD(4) * t241 + qJD(1);
	t255 = qJD(1) * t241 + qJD(4);
	t253 = t256 * t243;
	t252 = t241 * t271 + t269;
	t251 = -qJ(3) * t241 - t261 * t244 - pkin(1);
	t248 = t275 * t244;
	t247 = qJD(2) * t275;
	t246 = (t272 * t240 + t273 * t243) * qJD(4) - t276;
	t235 = t245 * t253 + (-t255 * t242 + t257) * t240;
	t234 = -t243 * t257 + t252 * qJD(4) + (t241 * t269 + t271) * qJD(1);
	t233 = t242 * t253 + (t255 * t245 + t258) * t240;
	t232 = -qJD(4) * t268 - t277 * t243 + t256 * t270;
	t1 = [t245 * t262 - t273 * t233 - t272 * t232 + t278 * t242 + (-t274 * t242 + t251 * t245) * qJD(1), (-t245 * t247 + t261 * t267) * t241 + (t246 * t245 - t267 * t275) * t244, -t241 * t267 + t257, t252 * qJD(5) - t273 * t234 + t272 * t235, t234, 0; t242 * t262 + t273 * t235 + t272 * t234 - t278 * t245 + (t251 * t242 + t274 * t245) * qJD(1), (-t261 * t241 + t248) * t266 + (-t241 * t247 + t246 * t244) * t242, t277, -(-t241 * t270 + t268) * qJD(5) + t272 * t233 - t273 * t232, t232, 0; 0, qJD(2) * t248 + t246 * t241, t265, (-t272 * t263 + t273 * t265) * t243 + (t272 * t265 + (t273 * qJD(4) - qJD(5)) * t244) * t240, -t240 * t263 - t243 * t265, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:02:42
	% EndTime: 2019-10-10 10:02:43
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (229->55), mult. (694->85), div. (0->0), fcn. (577->6), ass. (0->39)
	t205 = sin(qJ(2));
	t208 = cos(qJ(2));
	t207 = cos(qJ(4));
	t222 = pkin(2) + pkin(8) - r_i_i_C(3) - qJ(6);
	t240 = t222 * qJD(2) + qJD(5) * t207 - qJD(3);
	t245 = t240 * t205 - (pkin(3) + pkin(7)) * qJD(1) - (qJ(3) * qJD(2) - qJD(6)) * t208;
	t204 = sin(qJ(4));
	t226 = -r_i_i_C(1) - pkin(5) - pkin(4);
	t237 = r_i_i_C(2) + qJ(5);
	t215 = t226 * t204 + t237 * t207 - qJ(3);
	t244 = t222 * t205 + t215 * t208;
	t243 = t215 * qJD(2) + qJD(6);
	t206 = sin(qJ(1));
	t229 = qJD(2) * t208;
	t224 = t206 * t229;
	t209 = cos(qJ(1));
	t231 = qJD(1) * t209;
	t241 = t205 * t231 + t224;
	t236 = t206 * t204;
	t235 = t206 * t207;
	t234 = t209 * t204;
	t233 = t209 * t207;
	t232 = qJD(1) * t206;
	t230 = qJD(2) * t205;
	t228 = qJD(2) * t209;
	t227 = qJD(4) * t208;
	t223 = t208 * t228;
	t220 = qJD(4) * t205 + qJD(1);
	t219 = qJD(1) * t205 + qJD(4);
	t217 = t220 * t207;
	t216 = t205 * t234 + t235;
	t212 = t204 * qJD(5) + (-qJ(3) * t205 - t222 * t208 - pkin(1)) * qJD(1);
	t211 = (t237 * t204 - t226 * t207) * qJD(4) - t240;
	t210 = t243 * t205 + t211 * t208;
	t199 = t209 * t217 + (-t219 * t206 + t223) * t204;
	t198 = -t207 * t223 + t216 * qJD(4) + (t205 * t235 + t234) * qJD(1);
	t197 = t206 * t217 + (t219 * t209 + t224) * t204;
	t196 = -qJD(4) * t233 - t241 * t207 + t220 * t236;
	t1 = [-t237 * t196 + t226 * t197 + t245 * t206 + t212 * t209, t210 * t209 + t244 * t232, -t205 * t232 + t223, t216 * qJD(5) + t226 * t198 + t237 * t199, t198, t205 * t228 + t208 * t232; t237 * t198 - t226 * t199 + t212 * t206 - t245 * t209, t210 * t206 - t244 * t231, t241, -(-t205 * t236 + t233) * qJD(5) + t237 * t197 + t226 * t196, t196, t206 * t230 - t208 * t231; 0, t211 * t205 - t243 * t208, t230, (-t226 * t230 - t237 * t227) * t207 + (t237 * t230 + (-t226 * qJD(4) - qJD(5)) * t208) * t204, -t204 * t227 - t207 * t230, -t229;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end