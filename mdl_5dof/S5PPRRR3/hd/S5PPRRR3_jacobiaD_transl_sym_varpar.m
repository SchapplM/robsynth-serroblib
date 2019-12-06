% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PPRRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:17:26
	% EndTime: 2019-12-05 15:17:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:17:26
	% EndTime: 2019-12-05 15:17:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:17:26
	% EndTime: 2019-12-05 15:17:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:17:27
	% EndTime: 2019-12-05 15:17:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (30->18), div. (0->0), fcn. (24->6), ass. (0->10)
	t56 = sin(pkin(8));
	t59 = sin(qJ(3));
	t64 = t56 * t59;
	t60 = cos(qJ(3));
	t63 = t56 * t60;
	t58 = cos(pkin(8));
	t62 = t58 * t59;
	t61 = t58 * t60;
	t57 = cos(pkin(9));
	t1 = [0, 0, ((-t57 * t61 - t64) * r_i_i_C(1) + (t57 * t62 - t63) * r_i_i_C(2)) * qJD(3), 0, 0; 0, 0, ((-t57 * t63 + t62) * r_i_i_C(1) + (t57 * t64 + t61) * r_i_i_C(2)) * qJD(3), 0, 0; 0, 0, (-r_i_i_C(1) * t60 + r_i_i_C(2) * t59) * sin(pkin(9)) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:17:27
	% EndTime: 2019-12-05 15:17:27
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (55->26), mult. (191->56), div. (0->0), fcn. (172->8), ass. (0->28)
	t188 = cos(pkin(9));
	t189 = cos(pkin(8));
	t191 = sin(qJ(3));
	t202 = t189 * t191;
	t187 = sin(pkin(8));
	t193 = cos(qJ(3));
	t203 = t187 * t193;
	t209 = t188 * t202 - t203;
	t208 = -pkin(6) - r_i_i_C(3);
	t186 = sin(pkin(9));
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
	% StartTime: 2019-12-05 15:17:27
	% EndTime: 2019-12-05 15:17:27
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (180->39), mult. (333->73), div. (0->0), fcn. (304->10), ass. (0->40)
	t234 = cos(pkin(9));
	t235 = cos(pkin(8));
	t237 = sin(qJ(3));
	t257 = t235 * t237;
	t233 = sin(pkin(8));
	t239 = cos(qJ(3));
	t258 = t233 * t239;
	t267 = t234 * t257 - t258;
	t231 = qJ(4) + qJ(5);
	t228 = sin(t231);
	t229 = cos(t231);
	t230 = qJD(4) + qJD(5);
	t236 = sin(qJ(4));
	t266 = pkin(4) * qJD(4) * t236 + (r_i_i_C(1) * t228 + r_i_i_C(2) * t229) * t230;
	t265 = -r_i_i_C(3) - pkin(7) - pkin(6);
	t264 = t228 * t230;
	t263 = t229 * t230;
	t232 = sin(pkin(9));
	t262 = t230 * t232;
	t261 = t232 * t236;
	t260 = t232 * t239;
	t259 = t233 * t237;
	t256 = t235 * t239;
	t223 = t234 * t258 - t257;
	t244 = t234 * t259 + t256;
	t218 = t244 * qJD(3);
	t247 = -t233 * t262 + t218;
	t255 = (-t223 * t263 + t247 * t228) * r_i_i_C(1) + (t223 * t264 + t247 * t229) * r_i_i_C(2);
	t225 = t234 * t256 + t259;
	t220 = t267 * qJD(3);
	t246 = -t235 * t262 + t220;
	t254 = (-t225 * t263 + t246 * t228) * r_i_i_C(1) + (t225 * t264 + t246 * t229) * r_i_i_C(2);
	t249 = qJD(3) * t232 * t237;
	t243 = t230 * t234 + t249;
	t251 = t230 * t260;
	t253 = (t243 * t228 - t229 * t251) * r_i_i_C(1) + (t228 * t251 + t243 * t229) * r_i_i_C(2);
	t238 = cos(qJ(4));
	t245 = t238 * pkin(4) + r_i_i_C(1) * t229 - r_i_i_C(2) * t228 + pkin(3);
	t242 = qJD(3) * t245;
	t1 = [0, 0, t265 * t220 - t225 * t242 + t266 * t267, (t220 * t236 + (-t225 * t238 - t235 * t261) * qJD(4)) * pkin(4) + t254, t254; 0, 0, t265 * t218 - t223 * t242 + t244 * t266, (t218 * t236 + (-t223 * t238 - t233 * t261) * qJD(4)) * pkin(4) + t255, t255; 0, 0, (t266 * t237 + (t265 * t237 - t245 * t239) * qJD(3)) * t232, (t236 * t249 + (t234 * t236 - t238 * t260) * qJD(4)) * pkin(4) + t253, t253;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end