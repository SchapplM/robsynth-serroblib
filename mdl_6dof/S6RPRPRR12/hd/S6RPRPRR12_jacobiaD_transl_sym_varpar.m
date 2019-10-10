% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:06
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR12_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR12_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR12_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:05:58
	% EndTime: 2019-10-10 01:05:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:05:58
	% EndTime: 2019-10-10 01:05:58
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
	% StartTime: 2019-10-10 01:05:58
	% EndTime: 2019-10-10 01:05:58
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
	% StartTime: 2019-10-10 01:05:58
	% EndTime: 2019-10-10 01:05:58
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
	% StartTime: 2019-10-10 01:05:59
	% EndTime: 2019-10-10 01:05:59
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (48->21), mult. (142->35), div. (0->0), fcn. (98->4), ass. (0->17)
	t138 = sin(qJ(3));
	t140 = cos(qJ(3));
	t150 = r_i_i_C(3) + qJ(4);
	t151 = pkin(3) - r_i_i_C(2);
	t155 = -(t150 * t138 + t151 * t140) * qJD(3) + t140 * qJD(4);
	t154 = t151 * qJD(3) - qJD(4);
	t152 = t151 * t138 - t150 * t140 + qJ(2);
	t139 = sin(qJ(1));
	t149 = qJD(1) * t139;
	t141 = cos(qJ(1));
	t148 = qJD(1) * t141;
	t147 = qJD(3) * t139;
	t146 = qJD(3) * t141;
	t144 = qJD(1) * t151;
	t143 = qJD(1) * t150;
	t142 = qJD(2) + (-pkin(1) - pkin(7) - r_i_i_C(1)) * qJD(1) - t155;
	t1 = [t142 * t141 - t152 * t149, t148, (t141 * t144 + t150 * t147) * t140 + (-t154 * t139 + t141 * t143) * t138, t138 * t147 - t140 * t148, 0, 0; t142 * t139 + t152 * t148, t149, (t139 * t144 - t150 * t146) * t140 + (t139 * t143 + t154 * t141) * t138, -t138 * t146 - t140 * t149, 0, 0; 0, 0, t155, qJD(3) * t140, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:05:59
	% EndTime: 2019-10-10 01:05:59
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (105->44), mult. (326->75), div. (0->0), fcn. (254->6), ass. (0->34)
	t201 = cos(qJ(3));
	t216 = pkin(3) + pkin(8) + r_i_i_C(3);
	t227 = t216 * t201;
	t198 = sin(qJ(3));
	t226 = -qJ(4) * t201 + t216 * t198 + qJ(2);
	t197 = sin(qJ(5));
	t200 = cos(qJ(5));
	t205 = qJD(4) + (r_i_i_C(1) * t200 - r_i_i_C(2) * t197) * qJD(5);
	t225 = -t216 * qJD(3) + t205;
	t202 = cos(qJ(1));
	t224 = t200 * t202;
	t199 = sin(qJ(1));
	t223 = qJD(1) * t199;
	t222 = qJD(1) * t201;
	t221 = qJD(1) * t202;
	t220 = qJD(3) * t199;
	t219 = qJD(3) * t201;
	t218 = qJD(3) * t202;
	t217 = qJD(5) * t198;
	t215 = t198 * t220;
	t214 = t198 * t218;
	t212 = qJD(1) * t216;
	t211 = qJD(5) * t201 + qJD(1);
	t210 = qJD(5) + t222;
	t208 = t211 * t197;
	t207 = r_i_i_C(1) * t197 + r_i_i_C(2) * t200 + qJ(4);
	t206 = qJD(1) * t207;
	t204 = t210 * t199 + t214;
	t203 = -t201 * qJD(4) + qJD(2) + (-pkin(1) - pkin(4) - pkin(7)) * qJD(1) + (qJ(4) * t198 + t227) * qJD(3);
	t196 = t204 * t197 - t211 * t224;
	t195 = t204 * t200 + t202 * t208;
	t194 = t211 * t200 * t199 + (t210 * t202 - t215) * t197;
	t193 = -t210 * t224 + (qJD(3) * t198 * t200 + t208) * t199;
	t1 = [t196 * r_i_i_C(1) + t195 * r_i_i_C(2) + t203 * t202 - t226 * t223, t221, (t202 * t212 + t207 * t220) * t201 + (t225 * t199 + t202 * t206) * t198, -t201 * t221 + t215, t193 * r_i_i_C(1) + t194 * r_i_i_C(2), 0; -t194 * r_i_i_C(1) + t193 * r_i_i_C(2) + t203 * t199 + t226 * t221, t223, (t199 * t212 - t207 * t218) * t201 + (t199 * t206 - t225 * t202) * t198, -t199 * t222 - t214, -t195 * r_i_i_C(1) + t196 * r_i_i_C(2), 0; 0, 0, t205 * t201 + (-t207 * t198 - t227) * qJD(3), t219, (-t197 * t219 - t200 * t217) * r_i_i_C(2) + (-t197 * t217 + t200 * t219) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:05:59
	% EndTime: 2019-10-10 01:05:59
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (283->61), mult. (498->92), div. (0->0), fcn. (390->8), ass. (0->52)
	t232 = qJD(5) + qJD(6);
	t237 = cos(qJ(5));
	t268 = pkin(5) * qJD(5);
	t249 = t237 * t268 + qJD(4);
	t233 = qJ(5) + qJ(6);
	t230 = sin(t233);
	t269 = r_i_i_C(2) * t230;
	t231 = cos(t233);
	t270 = r_i_i_C(1) * t231;
	t242 = (-t269 + t270) * t232 + t249;
	t259 = pkin(3) + r_i_i_C(3) + pkin(9) + pkin(8);
	t280 = -t259 * qJD(3) + t242;
	t235 = sin(qJ(3));
	t238 = cos(qJ(3));
	t234 = sin(qJ(5));
	t255 = pkin(5) * t234 + qJ(4);
	t258 = t234 * t268;
	t279 = -t258 + (t259 * t235 - t255 * t238 + qJ(2)) * qJD(1);
	t239 = cos(qJ(1));
	t254 = t232 * t238 + qJD(1);
	t278 = t239 * t254;
	t277 = (qJD(5) * t238 + qJD(1)) * t234;
	t276 = t259 * t238;
	t264 = qJD(1) * t238;
	t253 = t232 + t264;
	t236 = sin(qJ(1));
	t262 = qJD(3) * t236;
	t256 = t235 * t262;
	t274 = t253 * t239 - t256;
	t272 = pkin(5) * t237;
	t271 = r_i_i_C(1) * t230;
	t267 = t232 * t235;
	t246 = t254 * t236;
	t224 = t230 * t246 - t274 * t231;
	t225 = t274 * t230 + t231 * t246;
	t266 = t224 * r_i_i_C(1) + t225 * r_i_i_C(2);
	t260 = qJD(3) * t239;
	t257 = t235 * t260;
	t243 = t253 * t236 + t257;
	t226 = t230 * t278 + t243 * t231;
	t227 = t243 * t230 - t231 * t278;
	t265 = -t226 * r_i_i_C(1) + t227 * r_i_i_C(2);
	t263 = qJD(1) * t239;
	t261 = qJD(3) * t238;
	t251 = -qJD(5) - t264;
	t250 = qJD(1) * t259;
	t247 = -r_i_i_C(2) * t231 - t271;
	t245 = -t247 + t255;
	t244 = qJD(1) * t245;
	t241 = qJD(2) - t249 * t238 + (-pkin(1) - pkin(7) - pkin(4) - t272) * qJD(1) + (t255 * t235 + t276) * qJD(3);
	t228 = t261 * t270;
	t1 = [t227 * r_i_i_C(1) + t226 * r_i_i_C(2) - t279 * t236 + t241 * t239, t263, (t239 * t250 + t245 * t262) * t238 + (t280 * t236 + t239 * t244) * t235, -t238 * t263 + t256, (t251 * t239 * t237 + (qJD(3) * t235 * t237 + t277) * t236) * pkin(5) + t266, t266; -t225 * r_i_i_C(1) + t224 * r_i_i_C(2) + t241 * t236 + t279 * t239, qJD(1) * t236, (t236 * t250 - t245 * t260) * t238 + (t236 * t244 - t280 * t239) * t235, -t236 * t264 - t257, (-t239 * t277 + (t251 * t236 - t257) * t237) * pkin(5) + t265, t265; 0, 0, t242 * t238 + (-t245 * t235 - t276) * qJD(3), t261, t228 + (-t269 + t272) * t261 + (t247 * t232 - t258) * t235, -t267 * t271 + t228 + (-t230 * t261 - t231 * t267) * r_i_i_C(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end