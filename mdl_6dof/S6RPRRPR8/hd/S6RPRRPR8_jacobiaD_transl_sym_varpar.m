% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:37
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
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:37
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:38
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 01:35:39
	% EndTime: 2019-10-10 01:35:39
	% DurationCPUTime: 0.23s
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
	% StartTime: 2019-10-10 01:35:39
	% EndTime: 2019-10-10 01:35:39
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (193->50), mult. (378->76), div. (0->0), fcn. (296->8), ass. (0->43)
	t213 = cos(qJ(4));
	t245 = t213 * pkin(4);
	t205 = pkin(3) + t245;
	t211 = sin(qJ(3));
	t214 = cos(qJ(3));
	t244 = r_i_i_C(3) + qJ(5) + pkin(8);
	t249 = t244 * t214;
	t256 = (-t205 * t211 - qJ(2) + t249) * qJD(1) - qJD(4) * t245;
	t208 = qJ(4) + pkin(10);
	t206 = sin(t208);
	t207 = cos(t208);
	t225 = r_i_i_C(1) * t207 - r_i_i_C(2) * t206;
	t222 = t205 + t225;
	t255 = -(-t222 * t211 + t249) * qJD(3) - qJD(5) * t211;
	t230 = t244 * t211;
	t254 = t222 * t214 + t230;
	t215 = cos(qJ(1));
	t238 = qJD(4) * t211;
	t227 = qJD(1) + t238;
	t223 = t227 * t215;
	t212 = sin(qJ(1));
	t226 = qJD(1) * t211 + qJD(4);
	t239 = qJD(3) * t215;
	t248 = t226 * t212 - t214 * t239;
	t240 = qJD(3) * t214;
	t247 = t212 * t240 + t226 * t215;
	t210 = sin(qJ(4));
	t246 = pkin(4) * t210;
	t243 = qJD(1) * t212;
	t242 = qJD(1) * t215;
	t241 = qJD(3) * t211;
	t237 = qJD(4) * t214;
	t235 = t214 * qJD(5);
	t224 = t227 * t212;
	t221 = r_i_i_C(1) * t206 + r_i_i_C(2) * t207 + t246;
	t220 = qJD(4) * t221;
	t217 = qJD(1) * t254;
	t216 = -t238 * t246 - t235 + qJD(2) + (t205 * t214 + t230) * qJD(3) + (-pkin(1) - pkin(7) - t246) * qJD(1);
	t204 = -t206 * t224 + t247 * t207;
	t203 = -t206 * t247 - t207 * t224;
	t202 = -t206 * t223 - t207 * t248;
	t201 = t248 * t206 - t207 * t223;
	t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t256 * t212 + t216 * t215, t242, t215 * t217 + (-t221 * t237 - t255) * t212, t203 * r_i_i_C(1) - t204 * r_i_i_C(2) + (-t210 * t247 - t213 * t224) * pkin(4), t212 * t241 - t214 * t242, 0; t204 * r_i_i_C(1) + t203 * r_i_i_C(2) + t216 * t212 - t256 * t215, t243, t212 * t217 + (t214 * t220 + t255) * t215, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2) + (-t210 * t248 + t213 * t223) * pkin(4), -t211 * t239 - t214 * t243, 0; 0, 0, -t254 * qJD(3) + t211 * t220 + t235, (-t225 - t245) * t237 + t221 * t241, t240, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:39
	% EndTime: 2019-10-10 01:35:39
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (425->58), mult. (488->84), div. (0->0), fcn. (386->10), ass. (0->51)
	t246 = qJ(4) + pkin(10);
	t236 = pkin(5) * cos(t246) + cos(qJ(4)) * pkin(4);
	t230 = t236 * qJD(4);
	t232 = pkin(3) + t236;
	t248 = sin(qJ(3));
	t251 = cos(qJ(3));
	t282 = r_i_i_C(3) + pkin(9) + qJ(5) + pkin(8);
	t289 = t282 * t251;
	t299 = (-t232 * t248 - qJ(2) + t289) * qJD(1) - t230;
	t235 = sin(qJ(4)) * pkin(4) + pkin(5) * sin(t246);
	t229 = t235 * qJD(4);
	t245 = qJD(4) + qJD(6);
	t242 = qJ(6) + t246;
	t238 = sin(t242);
	t239 = cos(t242);
	t296 = r_i_i_C(1) * t238 + r_i_i_C(2) * t239;
	t256 = t296 * t245 + t229;
	t284 = r_i_i_C(2) * t238;
	t285 = r_i_i_C(1) * t239;
	t257 = t232 - t284 + t285;
	t298 = (-t257 * t248 + t289) * qJD(3) + qJD(5) * t248 - t256 * t251;
	t266 = t282 * t248;
	t297 = t257 * t251 + t266;
	t249 = sin(qJ(1));
	t279 = qJD(1) * t248;
	t262 = t245 + t279;
	t252 = cos(qJ(1));
	t274 = qJD(3) * t252;
	t288 = t262 * t249 - t251 * t274;
	t275 = qJD(3) * t251;
	t287 = t249 * t275 + t262 * t252;
	t263 = -t245 * t248 - qJD(1);
	t258 = t263 * t252;
	t225 = t288 * t238 + t239 * t258;
	t226 = t238 * t258 - t288 * t239;
	t281 = -t225 * r_i_i_C(1) + t226 * r_i_i_C(2);
	t259 = t263 * t249;
	t227 = -t287 * t238 + t239 * t259;
	t228 = t238 * t259 + t287 * t239;
	t280 = t227 * r_i_i_C(1) - t228 * r_i_i_C(2);
	t278 = qJD(1) * t249;
	t277 = qJD(1) * t252;
	t276 = qJD(3) * t248;
	t272 = t251 * qJD(5);
	t271 = t245 * t285;
	t270 = t251 * t245 * t284 + t296 * t276;
	t260 = -t235 * t279 - t229;
	t255 = qJD(1) * t236 + t230 * t248 + t235 * t275;
	t254 = qJD(1) * t297;
	t253 = -t272 - t248 * t229 + qJD(2) + (t232 * t251 + t266) * qJD(3) + (-pkin(1) - pkin(7) - t235) * qJD(1);
	t1 = [t226 * r_i_i_C(1) + t225 * r_i_i_C(2) + t299 * t249 + t253 * t252, t277, t298 * t249 + t252 * t254, -t255 * t249 + t260 * t252 + t280, t249 * t276 - t251 * t277, t280; t228 * r_i_i_C(1) + t227 * r_i_i_C(2) + t253 * t249 - t299 * t252, t278, t249 * t254 - t298 * t252, t260 * t249 + t255 * t252 + t281, -t248 * t274 - t251 * t278, t281; 0, 0, -t297 * qJD(3) + t256 * t248 + t272, t235 * t276 + (-t230 - t271) * t251 + t270, t275, -t251 * t271 + t270;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end