% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRRP3
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PPRRP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRRP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (30->18), div. (0->0), fcn. (24->6), ass. (0->10)
	t56 = sin(pkin(7));
	t59 = sin(qJ(3));
	t64 = t56 * t59;
	t60 = cos(qJ(3));
	t63 = t56 * t60;
	t58 = cos(pkin(7));
	t62 = t58 * t59;
	t61 = t58 * t60;
	t57 = cos(pkin(8));
	t1 = [0, 0, ((-t57 * t61 - t64) * r_i_i_C(1) + (t57 * t62 - t63) * r_i_i_C(2)) * qJD(3), 0, 0; 0, 0, ((-t57 * t63 + t62) * r_i_i_C(1) + (t57 * t64 + t61) * r_i_i_C(2)) * qJD(3), 0, 0; 0, 0, (-r_i_i_C(1) * t60 + r_i_i_C(2) * t59) * sin(pkin(8)) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:39
	% EndTime: 2019-12-05 15:11:39
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (55->26), mult. (191->56), div. (0->0), fcn. (172->8), ass. (0->28)
	t188 = cos(pkin(8));
	t189 = cos(pkin(7));
	t191 = sin(qJ(3));
	t202 = t189 * t191;
	t187 = sin(pkin(7));
	t193 = cos(qJ(3));
	t203 = t187 * t193;
	t209 = t188 * t202 - t203;
	t208 = -pkin(6) - r_i_i_C(3);
	t186 = sin(pkin(8));
	t190 = sin(qJ(4));
	t207 = t186 * t190;
	t192 = cos(qJ(4));
	t206 = t186 * t192;
	t205 = t186 * t193;
	t204 = t187 * t191;
	t201 = t189 * t193;
	t199 = t190 * r_i_i_C(1) + t192 * r_i_i_C(2);
	t198 = t192 * r_i_i_C(1) - t190 * r_i_i_C(2) + pkin(3);
	t197 = t199 * t191;
	t184 = t188 * t201 + t204;
	t182 = t188 * t203 - t202;
	t196 = t188 * t204 + t201;
	t195 = qJD(4) * t199;
	t194 = qJD(3) * t198;
	t179 = t209 * qJD(3);
	t177 = t196 * qJD(3);
	t1 = [0, 0, t208 * t179 - t184 * t194 + t209 * t195, t199 * t179 + ((-t184 * t192 - t189 * t207) * r_i_i_C(1) + (t184 * t190 - t189 * t206) * r_i_i_C(2)) * qJD(4), 0; 0, 0, t208 * t177 - t182 * t194 + t196 * t195, t199 * t177 + ((-t182 * t192 - t187 * t207) * r_i_i_C(1) + (t182 * t190 - t187 * t206) * r_i_i_C(2)) * qJD(4), 0; 0, 0, (qJD(4) * t197 + (t208 * t191 - t198 * t193) * qJD(3)) * t186, t186 * qJD(3) * t197 + ((t188 * t190 - t192 * t205) * r_i_i_C(1) + (t188 * t192 + t190 * t205) * r_i_i_C(2)) * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:39
	% EndTime: 2019-12-05 15:11:39
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (124->36), mult. (404->65), div. (0->0), fcn. (384->8), ass. (0->35)
	t242 = sin(qJ(4));
	t244 = cos(qJ(4));
	t262 = r_i_i_C(3) + qJ(5);
	t263 = r_i_i_C(1) + pkin(4);
	t268 = (t263 * t242 - t262 * t244) * qJD(4) - qJD(5) * t242;
	t240 = cos(pkin(8));
	t241 = cos(pkin(7));
	t243 = sin(qJ(3));
	t256 = t241 * t243;
	t239 = sin(pkin(7));
	t245 = cos(qJ(3));
	t257 = t239 * t245;
	t267 = t240 * t256 - t257;
	t265 = t262 * t242 + t263 * t244 + pkin(3);
	t264 = -pkin(6) - r_i_i_C(2);
	t238 = sin(pkin(8));
	t261 = t238 * t242;
	t260 = t238 * t244;
	t259 = t238 * t245;
	t258 = t239 * t243;
	t255 = t241 * t245;
	t252 = qJD(3) * t238 * t243;
	t234 = t240 * t257 - t256;
	t251 = t234 * t244 + t239 * t261;
	t236 = t240 * t255 + t258;
	t250 = t236 * t244 + t241 * t261;
	t249 = t240 * t242 - t244 * t259;
	t248 = t240 * t258 + t255;
	t247 = qJD(3) * t265;
	t231 = t267 * qJD(3);
	t229 = t248 * qJD(3);
	t227 = t249 * qJD(4) + t242 * t252;
	t225 = t250 * qJD(4) - t231 * t242;
	t223 = t251 * qJD(4) - t229 * t242;
	t1 = [0, 0, t264 * t231 - t236 * t247 + t267 * t268, t250 * qJD(5) + t262 * (-t231 * t244 + (-t236 * t242 + t241 * t260) * qJD(4)) - t263 * t225, t225; 0, 0, t264 * t229 - t234 * t247 + t248 * t268, t251 * qJD(5) + t262 * (-t229 * t244 + (-t234 * t242 + t239 * t260) * qJD(4)) - t263 * t223, t223; 0, 0, (t268 * t243 + (t264 * t243 - t265 * t245) * qJD(3)) * t238, -t249 * qJD(5) - t262 * (t244 * t252 + (t240 * t244 + t242 * t259) * qJD(4)) + t263 * t227, -t227;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end