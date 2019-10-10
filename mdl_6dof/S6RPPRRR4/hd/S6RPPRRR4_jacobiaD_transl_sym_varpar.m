% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:06
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
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
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) - r_i_i_C(1);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:32
	% EndTime: 2019-10-10 00:06:32
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (14->11), mult. (36->20), div. (0->0), fcn. (26->4), ass. (0->8)
	t58 = -pkin(1) - pkin(2);
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t55 = cos(pkin(10));
	t54 = sin(pkin(10));
	t53 = (t54 * t57 - t55 * t56) * qJD(1);
	t52 = (t54 * t56 + t55 * t57) * qJD(1);
	t1 = [-t52 * r_i_i_C(1) + t53 * r_i_i_C(2) + t57 * qJD(2) + (-qJ(2) * t56 + t58 * t57) * qJD(1), qJD(1) * t57, 0, 0, 0, 0; t53 * r_i_i_C(1) + t52 * r_i_i_C(2) + t56 * qJD(2) + (qJ(2) * t57 + t58 * t56) * qJD(1), qJD(1) * t56, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:32
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (49->23), mult. (138->41), div. (0->0), fcn. (116->6), ass. (0->17)
	t47 = sin(qJ(4));
	t49 = cos(qJ(4));
	t51 = (r_i_i_C(1) * t47 + r_i_i_C(2) * t49) * qJD(4);
	t57 = -pkin(1) - pkin(2);
	t56 = -pkin(7) - r_i_i_C(3);
	t55 = qJD(4) * t47;
	t54 = qJD(4) * t49;
	t45 = sin(pkin(10));
	t46 = cos(pkin(10));
	t48 = sin(qJ(1));
	t50 = cos(qJ(1));
	t42 = t45 * t50 - t48 * t46;
	t43 = t48 * t45 + t46 * t50;
	t52 = r_i_i_C(1) * t49 - r_i_i_C(2) * t47 + pkin(3);
	t41 = t42 * qJD(1);
	t40 = t43 * qJD(1);
	t1 = [t50 * qJD(2) + t56 * t41 - t42 * t51 - t52 * t40 + (-qJ(2) * t48 + t57 * t50) * qJD(1), qJD(1) * t50, 0, (-t41 * t49 + t43 * t55) * r_i_i_C(2) + (-t41 * t47 - t43 * t54) * r_i_i_C(1), 0, 0; t48 * qJD(2) + t56 * t40 - t43 * t51 + t52 * t41 + (qJ(2) * t50 + t57 * t48) * qJD(1), qJD(1) * t48, 0, (-t40 * t49 - t42 * t55) * r_i_i_C(2) + (-t40 * t47 + t42 * t54) * r_i_i_C(1), 0, 0; 0, 0, 0, t51, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:33
	% EndTime: 2019-10-10 00:06:33
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (171->48), mult. (482->79), div. (0->0), fcn. (461->8), ass. (0->37)
	t235 = sin(qJ(4));
	t234 = sin(qJ(5));
	t236 = cos(qJ(5));
	t243 = t236 * r_i_i_C(1) - t234 * r_i_i_C(2) + pkin(4);
	t237 = cos(qJ(4));
	t260 = pkin(8) + r_i_i_C(3);
	t250 = t260 * t237;
	t239 = -t243 * t235 + t250;
	t267 = t239 * qJD(4);
	t249 = t260 * t235;
	t266 = t243 * t237 + t249;
	t256 = sin(pkin(10));
	t257 = cos(pkin(10));
	t258 = sin(qJ(1));
	t259 = cos(qJ(1));
	t229 = t259 * t256 - t258 * t257;
	t263 = qJD(1) * t258;
	t262 = qJD(1) * t259;
	t261 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
	t255 = qJD(4) * t235;
	t254 = qJD(4) * t237;
	t253 = qJD(5) * t234;
	t252 = qJD(5) * t236;
	t251 = qJD(5) * t237;
	t246 = t234 * r_i_i_C(1) + t236 * r_i_i_C(2);
	t228 = -t258 * t256 - t259 * t257;
	t226 = t228 * qJD(1);
	t245 = t228 * t251 + t226;
	t227 = t229 * qJD(1);
	t244 = t229 * t251 + t227;
	t242 = qJD(5) * t246;
	t241 = qJD(5) * t228 + t226 * t237 - t229 * t255;
	t240 = qJD(5) * t229 + t227 * t237 + t228 * t255;
	t238 = t266 * qJD(4) - t235 * t242;
	t225 = t245 * t234 + t240 * t236;
	t224 = -t240 * t234 + t245 * t236;
	t1 = [(-t227 * t234 + t228 * t252) * r_i_i_C(1) + (-t227 * t236 - t228 * t253) * r_i_i_C(2) - t227 * pkin(7) - qJ(2) * t263 - (-pkin(3) - t266) * t226 + (-t237 * t242 + t267) * t229 + t261 * t259, t262, 0, t239 * t227 + t238 * t228, t224 * r_i_i_C(1) - t225 * r_i_i_C(2), 0; t226 * pkin(7) + t225 * r_i_i_C(1) + t224 * r_i_i_C(2) + (pkin(4) * t235 - t250) * t228 * qJD(4) + qJ(2) * t262 + (pkin(4) * t237 + pkin(3) + t249) * t227 + t261 * t258, t263, 0, -t239 * t226 + t238 * t229, (t244 * r_i_i_C(1) + t241 * r_i_i_C(2)) * t236 + (t241 * r_i_i_C(1) - t244 * r_i_i_C(2)) * t234, 0; 0, 0, 0, t246 * t251 - t267, (-t235 * t253 + t236 * t254) * r_i_i_C(2) + (t234 * t254 + t235 * t252) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:33
	% EndTime: 2019-10-10 00:06:33
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (387->65), mult. (722->97), div. (0->0), fcn. (696->10), ass. (0->50)
	t280 = cos(qJ(4));
	t278 = sin(qJ(4));
	t279 = cos(qJ(5));
	t272 = pkin(5) * t279 + pkin(4);
	t276 = qJ(5) + qJ(6);
	t273 = sin(t276);
	t274 = cos(t276);
	t323 = r_i_i_C(1) * t274 - r_i_i_C(2) * t273;
	t290 = t272 + t323;
	t310 = r_i_i_C(3) + pkin(9) + pkin(8);
	t297 = t310 * t280;
	t283 = -t290 * t278 + t297;
	t275 = qJD(5) + qJD(6);
	t277 = sin(qJ(5));
	t298 = pkin(5) * qJD(5) * t277;
	t324 = t273 * r_i_i_C(1) + t274 * r_i_i_C(2);
	t318 = t275 * t324 + t298;
	t326 = qJD(4) * t283 - t318 * t280;
	t296 = t310 * t278;
	t325 = t290 * t280 + t296;
	t308 = sin(pkin(10));
	t309 = cos(pkin(10));
	t315 = sin(qJ(1));
	t316 = cos(qJ(1));
	t264 = t316 * t308 - t315 * t309;
	t320 = qJD(1) * t315;
	t319 = qJD(1) * t316;
	t317 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
	t263 = -t308 * t315 - t309 * t316;
	t307 = t263 * t275;
	t305 = t275 * t280;
	t261 = t263 * qJD(1);
	t302 = qJD(4) * t278;
	t289 = t261 * t280 - t264 * t302;
	t286 = t289 + t307;
	t262 = t264 * qJD(1);
	t291 = t264 * t305 + t262;
	t304 = (t273 * t286 + t274 * t291) * r_i_i_C(1) + (-t273 * t291 + t274 * t286) * r_i_i_C(2);
	t288 = -t262 * t280 - t263 * t302;
	t285 = t264 * t275 - t288;
	t292 = t263 * t305 + t261;
	t259 = -t273 * t285 + t274 * t292;
	t260 = t273 * t292 + t274 * t285;
	t303 = t259 * r_i_i_C(1) - t260 * r_i_i_C(2);
	t301 = qJD(4) * t280;
	t300 = qJD(5) * t279;
	t299 = qJD(5) * t280;
	t287 = t275 * t278 * t323 + t301 * t324;
	t282 = qJD(4) * t325 - t278 * t318;
	t1 = [(-t262 * t273 + t274 * t307) * r_i_i_C(1) + (-t262 * t274 - t273 * t307) * r_i_i_C(2) - t262 * pkin(7) + (-t262 * t277 + t263 * t300) * pkin(5) - qJ(2) * t320 - (-pkin(3) - t325) * t261 + t326 * t264 + t317 * t316, t319, 0, t262 * t283 + t263 * t282, ((t263 * t299 + t261) * t279 + (-qJD(5) * t264 + t288) * t277) * pkin(5) + t303, t303; t261 * pkin(7) + t260 * r_i_i_C(1) + t259 * r_i_i_C(2) + (t261 * t277 + t264 * t300) * pkin(5) + qJ(2) * t319 + (t280 * t298 + (t272 * t278 - t297) * qJD(4)) * t263 + (t272 * t280 + pkin(3) + t296) * t262 + t317 * t315, t320, 0, -t261 * t283 + t264 * t282, ((t264 * t299 + t262) * t279 + (qJD(5) * t263 + t289) * t277) * pkin(5) + t304, t304; 0, 0, 0, -t326, (t277 * t301 + t278 * t300) * pkin(5) + t287, t287;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end