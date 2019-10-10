% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:21
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
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
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) + r_i_i_C(2);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (23->16), mult. (72->29), div. (0->0), fcn. (46->4), ass. (0->13)
	t17 = sin(qJ(3));
	t19 = cos(qJ(3));
	t29 = (r_i_i_C(1) * t19 - r_i_i_C(2) * t17) * qJD(3);
	t18 = sin(qJ(1));
	t28 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t27 = qJD(1) * t20;
	t26 = qJD(3) * t18;
	t25 = qJD(3) * t20;
	t24 = -pkin(1) - pkin(7) - r_i_i_C(3);
	t22 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19 + qJ(2);
	t21 = qJD(2) + t29;
	t1 = [t21 * t20 + (-t18 * t22 + t20 * t24) * qJD(1), t27, (-t17 * t27 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 + t19 * t27) * r_i_i_C(1), 0, 0, 0; t21 * t18 + (t18 * t24 + t20 * t22) * qJD(1), t28, (-t17 * t28 + t19 * t25) * r_i_i_C(2) + (t17 * t25 + t19 * t28) * r_i_i_C(1), 0, 0, 0; 0, 0, -t29, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:50
	% EndTime: 2019-10-10 01:21:50
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (87->40), mult. (278->71), div. (0->0), fcn. (217->6), ass. (0->32)
	t201 = sin(qJ(3));
	t204 = cos(qJ(3));
	t223 = pkin(8) + r_i_i_C(3);
	t215 = t223 * t204;
	t229 = -pkin(3) * t201 - qJ(2) + t215;
	t200 = sin(qJ(4));
	t203 = cos(qJ(4));
	t209 = r_i_i_C(1) * t203 - r_i_i_C(2) * t200 + pkin(3);
	t216 = t223 * t201;
	t228 = t209 * t204 + t216;
	t205 = cos(qJ(1));
	t211 = qJD(1) * t201 + qJD(4);
	t226 = t211 * t205;
	t202 = sin(qJ(1));
	t218 = qJD(3) * t205;
	t225 = t211 * t202 - t204 * t218;
	t224 = -pkin(1) - pkin(7);
	t222 = qJD(1) * t202;
	t221 = qJD(1) * t205;
	t220 = qJD(3) * t201;
	t219 = qJD(3) * t204;
	t217 = qJD(4) * t204;
	t212 = -qJD(4) * t201 - qJD(1);
	t210 = r_i_i_C(1) * t200 + r_i_i_C(2) * t203;
	t208 = t212 * t205;
	t207 = t210 * qJD(4);
	t206 = qJD(2) + (pkin(3) * t204 + t216) * qJD(3);
	t199 = t203 * t226 + (t212 * t200 + t203 * t219) * t202;
	t198 = t212 * t203 * t202 + (-t202 * t219 - t226) * t200;
	t197 = t200 * t208 - t225 * t203;
	t196 = t225 * t200 + t203 * t208;
	t1 = [t197 * r_i_i_C(1) + t196 * r_i_i_C(2) + t206 * t205 + (t229 * t202 + t224 * t205) * qJD(1), t221, t228 * t221 + (-t210 * t217 + (-t209 * t201 + t215) * qJD(3)) * t202, t198 * r_i_i_C(1) - t199 * r_i_i_C(2), 0, 0; t199 * r_i_i_C(1) + t198 * r_i_i_C(2) + t206 * t202 + (t224 * t202 - t229 * t205) * qJD(1), t222, (t209 * t218 + t223 * t222) * t201 + (t209 * t222 + (-t223 * qJD(3) + t207) * t205) * t204, -t196 * r_i_i_C(1) + t197 * r_i_i_C(2), 0, 0; 0, 0, -t228 * qJD(3) + t201 * t207, (t200 * t217 + t203 * t220) * r_i_i_C(2) + (t200 * t220 - t203 * t217) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:50
	% EndTime: 2019-10-10 01:21:51
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (168->55), mult. (518->86), div. (0->0), fcn. (433->6), ass. (0->38)
	t241 = cos(qJ(3));
	t237 = sin(qJ(4));
	t240 = cos(qJ(4));
	t266 = r_i_i_C(3) + qJ(5);
	t267 = r_i_i_C(2) - pkin(4);
	t246 = t266 * t237 - t267 * t240 + pkin(3);
	t238 = sin(qJ(3));
	t268 = pkin(8) + r_i_i_C(1);
	t254 = t268 * t238;
	t273 = t246 * t241 + t254;
	t253 = t268 * t241;
	t272 = -pkin(3) * t238 - qJ(2) + t253;
	t256 = qJD(5) * t237;
	t271 = t256 + (t267 * t237 + t266 * t240) * qJD(4);
	t269 = -pkin(1) - pkin(7);
	t239 = sin(qJ(1));
	t265 = t239 * t237;
	t264 = t239 * t240;
	t242 = cos(qJ(1));
	t263 = t242 * t237;
	t262 = qJD(1) * t239;
	t261 = qJD(1) * t242;
	t260 = qJD(3) * t238;
	t259 = qJD(3) * t242;
	t258 = qJD(4) * t238;
	t257 = qJD(4) * t241;
	t255 = t240 * qJD(5);
	t252 = t242 * t238 * t240;
	t251 = t241 * t259;
	t249 = qJD(1) * t238 + qJD(4);
	t248 = t238 * t264 + t263;
	t245 = t239 * qJD(3) * t241 + t249 * t242;
	t243 = t238 * t256 + qJD(2) + (pkin(3) * t241 + t254) * qJD(3);
	t230 = -t237 * t262 + t245 * t240 - t258 * t265;
	t229 = (qJD(1) + t258) * t264 + t245 * t237;
	t228 = -t240 * t251 + (t238 * t263 + t264) * qJD(4) + t248 * qJD(1);
	t227 = -qJD(4) * t252 - t237 * t251 - t240 * t261 + t249 * t265;
	t1 = [t239 * t255 + t267 * t228 - t266 * t227 + t243 * t242 + (t272 * t239 + t269 * t242) * qJD(1), t261, t273 * t261 + (t271 * t241 + (-t238 * t246 + t253) * qJD(3)) * t239, t248 * qJD(5) + t267 * t229 + t266 * t230, t229, 0; -t242 * t255 - t267 * t230 + t266 * t229 + t243 * t239 + (t269 * t239 - t272 * t242) * qJD(1), t262, (t246 * t259 + t268 * t262) * t238 + (t246 * t262 + (-t268 * qJD(3) - t271) * t242) * t241, -(t252 - t265) * qJD(5) + t266 * t228 + t267 * t227, t227, 0; 0, 0, -t273 * qJD(3) - t238 * t271, (-t266 * t257 - t267 * t260) * t237 + (-t266 * t260 + (t267 * qJD(4) + qJD(5)) * t241) * t240, -t237 * t260 + t240 * t257, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:50
	% EndTime: 2019-10-10 01:21:51
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (228->59), mult. (692->91), div. (0->0), fcn. (588->6), ass. (0->41)
	t241 = sin(qJ(4));
	t244 = cos(qJ(4));
	t263 = pkin(4) + r_i_i_C(3) + qJ(6);
	t274 = r_i_i_C(2) + qJ(5);
	t252 = t263 * t241 - t274 * t244;
	t264 = pkin(5) + pkin(8) + r_i_i_C(1);
	t280 = t264 * qJD(3) - t252 * qJD(4) + qJD(5) * t241 + qJD(6) * t244;
	t253 = -t274 * t241 - t263 * t244;
	t279 = -pkin(3) + t253;
	t242 = sin(qJ(3));
	t245 = cos(qJ(3));
	t277 = -pkin(3) * t242 + t264 * t245 - qJ(2);
	t275 = -pkin(1) - pkin(7);
	t243 = sin(qJ(1));
	t273 = t243 * t241;
	t272 = t243 * t244;
	t246 = cos(qJ(1));
	t271 = t246 * t241;
	t270 = t246 * t244;
	t269 = qJD(1) * t243;
	t268 = qJD(1) * t246;
	t267 = qJD(3) * t242;
	t266 = qJD(3) * t246;
	t265 = qJD(4) * t245;
	t262 = t242 * t273;
	t261 = t242 * t270;
	t260 = t245 * t266;
	t258 = t264 * t242;
	t257 = qJD(1) * t242 + qJD(4);
	t255 = t242 * t272 + t271;
	t254 = t242 * t271 + t272;
	t250 = qJD(3) * t243 * t245 + t257 * t246;
	t249 = qJD(2) + (pkin(3) * t245 + t258) * qJD(3);
	t248 = qJD(3) * t279;
	t234 = t261 - t273;
	t231 = t262 - t270;
	t230 = -qJD(4) * t262 - t241 * t269 + t250 * t244;
	t229 = (qJD(4) * t242 + qJD(1)) * t272 + t250 * t241;
	t228 = t255 * qJD(1) + t254 * qJD(4) - t244 * t260;
	t227 = -qJD(4) * t261 - t241 * t260 - t244 * t268 + t257 * t273;
	t1 = [t254 * qJD(5) + t234 * qJD(6) - t274 * t227 - t263 * t228 + t249 * t246 + (t243 * t277 + t275 * t246) * qJD(1), t268, (-t245 * t279 + t258) * t268 + (t242 * t248 + t280 * t245) * t243, qJD(5) * t255 - t231 * qJD(6) - t263 * t229 + t274 * t230, t229, t230; t231 * qJD(5) + t255 * qJD(6) + t274 * t229 + t263 * t230 + t249 * t243 + (t275 * t243 - t246 * t277) * qJD(1), t269, (t264 * t269 - t266 * t279) * t242 + (-t246 * t280 - t269 * t279) * t245, -t234 * qJD(5) + qJD(6) * t254 - t263 * t227 + t274 * t228, t227, t228; 0, 0, -t242 * t280 + t245 * t248, t252 * t267 + (t253 * qJD(4) + qJD(5) * t244 - qJD(6) * t241) * t245, -t241 * t267 + t244 * t265, -t241 * t265 - t244 * t267;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end