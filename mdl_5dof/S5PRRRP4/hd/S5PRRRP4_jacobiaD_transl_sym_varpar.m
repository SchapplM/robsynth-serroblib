% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:33
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRP4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRP4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:33:46
	% EndTime: 2019-10-24 10:33:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:33:46
	% EndTime: 2019-10-24 10:33:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:33:46
	% EndTime: 2019-10-24 10:33:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(8)) * t14, 0, 0, 0; 0, sin(pkin(8)) * t14, 0, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:33:46
	% EndTime: 2019-10-24 10:33:46
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (33->9), mult. (40->14), div. (0->0), fcn. (25->6), ass. (0->14)
	t35 = pkin(2) * qJD(2);
	t27 = qJ(2) + qJ(3);
	t25 = cos(t27);
	t26 = qJD(2) + qJD(3);
	t34 = r_i_i_C(1) * t25 * t26;
	t24 = sin(t27);
	t33 = r_i_i_C(2) * t24 * t26;
	t32 = (-r_i_i_C(1) * t24 - r_i_i_C(2) * t25) * t26;
	t31 = -cos(qJ(2)) * t35 - t34;
	t29 = cos(pkin(8));
	t28 = sin(pkin(8));
	t23 = t29 * t33;
	t22 = t28 * t33;
	t1 = [0, t31 * t29 + t23, -t29 * t34 + t23, 0, 0; 0, t31 * t28 + t22, -t28 * t34 + t22, 0, 0; 0, -sin(qJ(2)) * t35 + t32, t32, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:33:46
	% EndTime: 2019-10-24 10:33:46
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (142->29), mult. (192->54), div. (0->0), fcn. (139->8), ass. (0->29)
	t189 = qJ(2) + qJ(3);
	t186 = sin(t189);
	t187 = cos(t189);
	t193 = cos(qJ(4));
	t207 = qJD(4) * t193;
	t188 = qJD(2) + qJD(3);
	t192 = sin(qJ(4));
	t214 = t188 * t192;
	t220 = t186 * t207 + t187 * t214;
	t219 = pkin(7) + r_i_i_C(3);
	t208 = qJD(4) * t192;
	t203 = t186 * t208;
	t218 = r_i_i_C(1) * t203 + t220 * r_i_i_C(2);
	t217 = pkin(2) * qJD(2);
	t216 = t186 * t188;
	t213 = t188 * t193;
	t190 = sin(pkin(8));
	t212 = t190 * t192;
	t211 = t190 * t193;
	t191 = cos(pkin(8));
	t210 = t191 * t192;
	t209 = t191 * t193;
	t205 = t218 * t190;
	t204 = t218 * t191;
	t198 = (r_i_i_C(1) * t192 + r_i_i_C(2) * t193) * t216;
	t197 = ((-r_i_i_C(1) * t193 - pkin(3)) * t187 - t219 * t186) * t188;
	t196 = -cos(qJ(2)) * t217 + t197;
	t195 = -pkin(3) * t216 + (-t186 * t213 - t187 * t208) * r_i_i_C(1) + t219 * t187 * t188 + (t186 * t214 - t187 * t207) * r_i_i_C(2);
	t1 = [0, t196 * t191 + t204, t191 * t197 + t204, t191 * t198 + ((-t187 * t209 - t212) * r_i_i_C(1) + (t187 * t210 - t211) * r_i_i_C(2)) * qJD(4), 0; 0, t196 * t190 + t205, t190 * t197 + t205, t190 * t198 + ((-t187 * t211 + t210) * r_i_i_C(1) + (t187 * t212 + t209) * r_i_i_C(2)) * qJD(4), 0; 0, -sin(qJ(2)) * t217 + t195, t195, (-t187 * t213 + t203) * r_i_i_C(2) - t220 * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:33:47
	% EndTime: 2019-10-24 10:33:47
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (270->37), mult. (374->58), div. (0->0), fcn. (291->8), ass. (0->38)
	t253 = sin(qJ(4));
	t283 = pkin(4) + r_i_i_C(1);
	t269 = t283 * t253;
	t287 = qJD(4) * t269;
	t254 = cos(qJ(4));
	t271 = qJD(4) * t254;
	t282 = r_i_i_C(3) + qJ(5);
	t286 = -qJD(5) * t253 - t282 * t271;
	t285 = pkin(7) + r_i_i_C(2);
	t250 = qJ(2) + qJ(3);
	t247 = sin(t250);
	t284 = t247 * t287;
	t281 = pkin(2) * qJD(2);
	t249 = qJD(2) + qJD(3);
	t280 = t247 * t249;
	t248 = cos(t250);
	t279 = t248 * t249;
	t278 = t249 * t253;
	t251 = sin(pkin(8));
	t277 = t251 * t253;
	t276 = t251 * t254;
	t252 = cos(pkin(8));
	t275 = t252 * t253;
	t274 = t252 * t254;
	t273 = t251 * t284;
	t272 = t252 * t284;
	t268 = t247 * t278;
	t267 = t254 * t280;
	t262 = -t248 * t274 - t277;
	t261 = -t248 * t276 + t275;
	t260 = -t282 * t253 - t283 * t254;
	t259 = -pkin(3) + t260;
	t258 = t259 * t280 + t285 * t279 + (-t286 - t287) * t248;
	t257 = t286 * t247 + (-t285 * t247 + t259 * t248) * t249;
	t256 = -cos(qJ(2)) * t281 + t257;
	t236 = t262 * qJD(4) + t252 * t268;
	t234 = t261 * qJD(4) + t251 * t268;
	t1 = [0, t256 * t252 + t272, t257 * t252 + t272, -t262 * qJD(5) - t282 * (t252 * t267 + (t248 * t275 - t276) * qJD(4)) + t283 * t236, -t236; 0, t256 * t251 + t273, t257 * t251 + t273, -t261 * qJD(5) - t282 * (t251 * t267 + (t248 * t277 + t274) * qJD(4)) + t283 * t234, -t234; 0, -sin(qJ(2)) * t281 + t258, t258, (t282 * t254 - t269) * t279 + (t260 * qJD(4) + qJD(5) * t254) * t247, t247 * t271 + t248 * t278;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end