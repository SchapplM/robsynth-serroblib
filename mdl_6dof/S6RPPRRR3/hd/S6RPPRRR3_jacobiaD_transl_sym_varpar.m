% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
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
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(10);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->9), mult. (24->12), div. (0->0), fcn. (14->4), ass. (0->6)
	t13 = -pkin(2) + r_i_i_C(2);
	t12 = r_i_i_C(3) + qJ(3);
	t11 = qJ(1) + pkin(10);
	t10 = cos(t11);
	t9 = sin(t11);
	t1 = [t10 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t12 * t9 + t13 * t10) * qJD(1), 0, qJD(1) * t10, 0, 0, 0; t9 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t13 * t9 + t12 * t10) * qJD(1), 0, qJD(1) * t9, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (53->19), mult. (76->33), div. (0->0), fcn. (48->6), ass. (0->14)
	t24 = sin(qJ(4));
	t25 = cos(qJ(4));
	t34 = (r_i_i_C(1) * t25 - r_i_i_C(2) * t24) * qJD(4);
	t33 = qJD(1) * t24;
	t32 = qJD(1) * t25;
	t31 = qJD(4) * t24;
	t30 = qJD(4) * t25;
	t29 = -pkin(2) - pkin(7) - r_i_i_C(3);
	t27 = r_i_i_C(1) * t24 + r_i_i_C(2) * t25 + qJ(3);
	t26 = qJD(3) + t34;
	t23 = qJ(1) + pkin(10);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t26 * t22 + (-cos(qJ(1)) * pkin(1) + t29 * t22 - t27 * t21) * qJD(1), 0, qJD(1) * t22, (-t21 * t30 - t22 * t33) * r_i_i_C(2) + (-t21 * t31 + t22 * t32) * r_i_i_C(1), 0, 0; t26 * t21 + (-sin(qJ(1)) * pkin(1) + t27 * t22 + t29 * t21) * qJD(1), 0, qJD(1) * t21, (-t21 * t33 + t22 * t30) * r_i_i_C(2) + (t21 * t32 + t22 * t31) * r_i_i_C(1), 0, 0; 0, 0, 0, -t34, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:50
	% EndTime: 2019-10-10 00:04:50
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (175->40), mult. (282->67), div. (0->0), fcn. (219->8), ass. (0->33)
	t202 = sin(qJ(4));
	t201 = sin(qJ(5));
	t203 = cos(qJ(5));
	t210 = r_i_i_C(1) * t203 - r_i_i_C(2) * t201 + pkin(4);
	t204 = cos(qJ(4));
	t221 = pkin(8) + r_i_i_C(3);
	t223 = t221 * t204;
	t231 = (-t210 * t202 + t223) * qJD(4);
	t216 = t221 * t202;
	t230 = t210 * t204 + t216;
	t229 = -pkin(4) * t202 - qJ(3) + t223;
	t212 = qJD(1) * t202 + qJD(5);
	t227 = t201 * t212;
	t226 = t203 * t212;
	t222 = -pkin(2) - pkin(7);
	t220 = qJD(4) * t202;
	t219 = qJD(4) * t204;
	t218 = qJD(5) * t204;
	t213 = -qJD(5) * t202 - qJD(1);
	t211 = r_i_i_C(1) * t201 + r_i_i_C(2) * t203;
	t209 = qJD(5) * t211;
	t208 = qJD(3) + (pkin(4) * t204 + t216) * qJD(4);
	t207 = -t201 * t219 + t213 * t203;
	t206 = t213 * t201 + t203 * t219;
	t205 = qJD(1) * t230;
	t200 = qJ(1) + pkin(10);
	t199 = cos(t200);
	t198 = sin(t200);
	t197 = t206 * t198 + t199 * t226;
	t196 = t207 * t198 - t199 * t227;
	t195 = -t198 * t226 + t206 * t199;
	t194 = t198 * t227 + t207 * t199;
	t1 = [t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t208 * t199 + (-cos(qJ(1)) * pkin(1) + t222 * t199 + t229 * t198) * qJD(1), 0, qJD(1) * t199, t199 * t205 + (-t211 * t218 + t231) * t198, t196 * r_i_i_C(1) - t197 * r_i_i_C(2), 0; t197 * r_i_i_C(1) + t196 * r_i_i_C(2) + t208 * t198 + (-sin(qJ(1)) * pkin(1) + t222 * t198 - t229 * t199) * qJD(1), 0, qJD(1) * t198, t198 * t205 + (t204 * t209 - t231) * t199, -t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0; 0, 0, 0, -t230 * qJD(4) + t202 * t209, (t201 * t218 + t203 * t220) * r_i_i_C(2) + (t201 * t220 - t203 * t218) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:50
	% EndTime: 2019-10-10 00:04:51
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (391->53), mult. (432->83), div. (0->0), fcn. (340->10), ass. (0->48)
	t238 = sin(qJ(4));
	t240 = cos(qJ(4));
	t234 = qJD(5) + qJD(6);
	t237 = sin(qJ(5));
	t274 = pkin(5) * t237;
	t261 = qJD(5) * t274;
	t236 = qJ(5) + qJ(6);
	t232 = sin(t236);
	t233 = cos(t236);
	t285 = r_i_i_C(1) * t232 + r_i_i_C(2) * t233;
	t244 = t285 * t234 + t261;
	t239 = cos(qJ(5));
	t229 = t239 * pkin(5) + pkin(4);
	t284 = -r_i_i_C(1) * t233 + r_i_i_C(2) * t232;
	t248 = t229 - t284;
	t269 = r_i_i_C(3) + pkin(9) + pkin(8);
	t276 = t269 * t240;
	t287 = (-t248 * t238 + t276) * qJD(4) - t244 * t240;
	t257 = t269 * t238;
	t286 = t248 * t240 + t257;
	t283 = -t229 * t238 - qJ(3) + t276;
	t265 = qJD(1) * t238;
	t253 = t234 + t265;
	t281 = t232 * t253;
	t280 = t233 * t253;
	t263 = qJD(4) * t240;
	t275 = (qJD(5) * t238 + qJD(1)) * t239 + t237 * t263;
	t235 = qJ(1) + pkin(10);
	t230 = sin(t235);
	t231 = cos(t235);
	t254 = -t234 * t238 - qJD(1);
	t246 = -t232 * t263 + t254 * t233;
	t222 = t230 * t281 + t246 * t231;
	t245 = t254 * t232 + t233 * t263;
	t223 = -t230 * t280 + t245 * t231;
	t267 = -t222 * r_i_i_C(1) + t223 * r_i_i_C(2);
	t224 = t246 * t230 - t231 * t281;
	t225 = t245 * t230 + t231 * t280;
	t266 = t224 * r_i_i_C(1) - t225 * r_i_i_C(2);
	t264 = qJD(4) * t238;
	t262 = qJD(5) * t239;
	t260 = pkin(5) * t262;
	t252 = -pkin(2) - pkin(7) - t274;
	t249 = (-qJD(5) - t265) * t237;
	t247 = t284 * t234 * t240 + t285 * t264;
	t243 = qJD(1) * t286;
	t242 = -t238 * t261 + qJD(3) + (t229 * t240 + t257) * qJD(4);
	t1 = [-t230 * t260 + t223 * r_i_i_C(1) + t222 * r_i_i_C(2) + t242 * t231 + (-cos(qJ(1)) * pkin(1) + t252 * t231 + t283 * t230) * qJD(1), 0, qJD(1) * t231, t287 * t230 + t231 * t243, (-t275 * t230 + t231 * t249) * pkin(5) + t266, t266; t231 * t260 + t225 * r_i_i_C(1) + t224 * r_i_i_C(2) + t242 * t230 + (-sin(qJ(1)) * pkin(1) + t252 * t230 - t283 * t231) * qJD(1), 0, qJD(1) * t230, t230 * t243 - t287 * t231, (t230 * t249 + t275 * t231) * pkin(5) + t267, t267; 0, 0, 0, -t286 * qJD(4) + t244 * t238, (t237 * t264 - t240 * t262) * pkin(5) + t247, t247;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end