% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:51
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
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
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
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
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (43->18), mult. (68->31), div. (0->0), fcn. (42->6), ass. (0->13)
	t24 = sin(qJ(3));
	t25 = cos(qJ(3));
	t26 = (r_i_i_C(1) * t24 + r_i_i_C(2) * t25) * qJD(3);
	t33 = pkin(7) + r_i_i_C(3);
	t32 = qJD(1) * t24;
	t31 = qJD(1) * t25;
	t30 = qJD(3) * t24;
	t29 = qJD(3) * t25;
	t27 = -r_i_i_C(1) * t25 + r_i_i_C(2) * t24 - pkin(2);
	t23 = qJ(1) + pkin(10);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t21 * t26 + (-cos(qJ(1)) * pkin(1) - t33 * t21 + t27 * t22) * qJD(1), 0, (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0, 0; -t22 * t26 + (-sin(qJ(1)) * pkin(1) + t33 * t22 + t27 * t21) * qJD(1), 0, (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t26, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:28
	% EndTime: 2019-10-10 00:51:28
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (92->20), mult. (138->30), div. (0->0), fcn. (94->6), ass. (0->17)
	t147 = sin(qJ(3));
	t148 = cos(qJ(3));
	t158 = r_i_i_C(3) + qJ(4);
	t160 = pkin(3) - r_i_i_C(2);
	t153 = t160 * t147 - t158 * t148;
	t162 = qJD(1) * t153;
	t161 = t153 * qJD(3) - t147 * qJD(4);
	t159 = pkin(7) + r_i_i_C(1);
	t157 = qJD(1) * t147;
	t156 = qJD(3) * t148;
	t154 = -t158 * t147 - t160 * t148;
	t151 = -pkin(2) + t154;
	t150 = t154 * qJD(3) + qJD(4) * t148;
	t146 = qJ(1) + pkin(10);
	t145 = cos(t146);
	t144 = sin(t146);
	t1 = [t161 * t144 + (-cos(qJ(1)) * pkin(1) - t159 * t144 + t151 * t145) * qJD(1), 0, t144 * t162 + t150 * t145, -t144 * t157 + t145 * t156, 0, 0; -t161 * t145 + (-sin(qJ(1)) * pkin(1) + t159 * t145 + t151 * t144) * qJD(1), 0, t150 * t144 - t145 * t162, t144 * t156 + t145 * t157, 0, 0; 0, 0, -t161, qJD(3) * t147, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:28
	% EndTime: 2019-10-10 00:51:28
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (201->42), mult. (322->69), div. (0->0), fcn. (250->8), ass. (0->33)
	t204 = cos(qJ(3));
	t201 = sin(qJ(5));
	t203 = cos(qJ(5));
	t209 = r_i_i_C(1) * t201 + r_i_i_C(2) * t203 + qJ(4);
	t202 = sin(qJ(3));
	t219 = pkin(3) + pkin(8) + r_i_i_C(3);
	t229 = t219 * t202;
	t231 = -t209 * t204 + t229;
	t234 = qJD(1) * t231;
	t233 = (-qJ(4) * t204 + t229) * qJD(3) - t202 * qJD(4);
	t213 = qJD(5) * t202 + qJD(1);
	t222 = qJD(3) * t204;
	t228 = t213 * t201 - t203 * t222;
	t227 = t201 * t222 + t213 * t203;
	t226 = pkin(4) + pkin(7);
	t224 = qJD(1) * t202;
	t223 = qJD(3) * t202;
	t221 = qJD(5) * t204;
	t215 = t219 * t204;
	t212 = -qJD(5) - t224;
	t211 = t212 * t201;
	t210 = t212 * t203;
	t208 = -qJ(4) * t202 - pkin(2) - t215;
	t207 = qJD(4) + (r_i_i_C(1) * t203 - r_i_i_C(2) * t201) * qJD(5);
	t205 = t207 * t204 + (-t209 * t202 - t215) * qJD(3);
	t200 = qJ(1) + pkin(10);
	t199 = cos(t200);
	t198 = sin(t200);
	t197 = t198 * t211 + t227 * t199;
	t196 = t198 * t210 - t228 * t199;
	t195 = -t227 * t198 + t199 * t211;
	t194 = t228 * t198 + t199 * t210;
	t1 = [t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t233 * t198 + (-cos(qJ(1)) * pkin(1) - t226 * t198 + t208 * t199) * qJD(1), 0, t198 * t234 + t205 * t199, -t198 * t224 + t199 * t222, t196 * r_i_i_C(1) - t197 * r_i_i_C(2), 0; t197 * r_i_i_C(1) + t196 * r_i_i_C(2) - t233 * t199 + (-sin(qJ(1)) * pkin(1) + t226 * t199 + t208 * t198) * qJD(1), 0, t205 * t198 - t199 * t234, t198 * t222 + t199 * t224, -t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0; 0, 0, -qJD(3) * t231 + t207 * t202, t223, (-t201 * t223 + t203 * t221) * r_i_i_C(2) + (t201 * t221 + t203 * t223) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:28
	% EndTime: 2019-10-10 00:51:29
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (423->55), mult. (494->85), div. (0->0), fcn. (386->10), ass. (0->48)
	t242 = cos(qJ(3));
	t239 = sin(qJ(5));
	t260 = pkin(5) * t239 + qJ(4);
	t238 = qJ(5) + qJ(6);
	t234 = sin(t238);
	t235 = cos(t238);
	t288 = r_i_i_C(1) * t234 + r_i_i_C(2) * t235;
	t249 = t260 + t288;
	t240 = sin(qJ(3));
	t264 = pkin(3) + r_i_i_C(3) + pkin(9) + pkin(8);
	t280 = t264 * t240;
	t285 = -t249 * t242 + t280;
	t290 = qJD(1) * t285;
	t241 = cos(qJ(5));
	t273 = t241 * pkin(5);
	t253 = qJD(5) * t273 + qJD(4);
	t289 = (-t260 * t242 + t280) * qJD(3) - t253 * t240;
	t287 = r_i_i_C(1) * t235 - r_i_i_C(2) * t234;
	t268 = qJD(1) * t240;
	t284 = t241 * (qJD(5) + t268);
	t236 = qJD(5) + qJD(6);
	t259 = t236 * t240 + qJD(1);
	t266 = qJD(3) * t242;
	t279 = t259 * t234 - t235 * t266;
	t278 = t234 * t266 + t259 * t235;
	t272 = pkin(7) + pkin(4) + t273;
	t237 = qJ(1) + pkin(10);
	t232 = sin(t237);
	t233 = cos(t237);
	t258 = -t236 - t268;
	t251 = t258 * t235;
	t224 = t279 * t232 + t233 * t251;
	t252 = t258 * t234;
	t225 = -t278 * t232 + t233 * t252;
	t270 = -t224 * r_i_i_C(1) + t225 * r_i_i_C(2);
	t226 = t232 * t251 - t279 * t233;
	t227 = t232 * t252 + t278 * t233;
	t269 = t226 * r_i_i_C(1) - t227 * r_i_i_C(2);
	t267 = qJD(3) * t240;
	t265 = qJD(5) * t239;
	t263 = pkin(5) * t265;
	t255 = t264 * t242;
	t250 = t288 * t236 * t242 + t287 * t267;
	t248 = t241 * t266 + (-qJD(5) * t240 - qJD(1)) * t239;
	t247 = -t260 * t240 - pkin(2) - t255;
	t246 = t287 * t236 + t253;
	t244 = t246 * t242 + (-t249 * t240 - t255) * qJD(3);
	t1 = [-t233 * t263 + t225 * r_i_i_C(1) + t224 * r_i_i_C(2) + t289 * t232 + (-cos(qJ(1)) * pkin(1) - t272 * t232 + t247 * t233) * qJD(1), 0, t232 * t290 + t244 * t233, -t232 * t268 + t233 * t266, (-t232 * t284 + t248 * t233) * pkin(5) + t269, t269; -t232 * t263 + t227 * r_i_i_C(1) + t226 * r_i_i_C(2) - t289 * t233 + (-sin(qJ(1)) * pkin(1) + t272 * t233 + t247 * t232) * qJD(1), 0, t244 * t232 - t233 * t290, t232 * t266 + t233 * t268, (t248 * t232 + t233 * t284) * pkin(5) + t270, t270; 0, 0, -qJD(3) * t285 + t246 * t240, t267, (t241 * t267 + t242 * t265) * pkin(5) + t250, t250;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end