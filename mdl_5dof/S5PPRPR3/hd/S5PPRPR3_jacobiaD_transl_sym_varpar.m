% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRPR3
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PPRPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
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
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->13), mult. (47->28), div. (0->0), fcn. (36->8), ass. (0->12)
	t64 = sin(pkin(7));
	t65 = cos(pkin(8));
	t71 = t64 * t65;
	t66 = cos(pkin(7));
	t70 = t65 * t66;
	t68 = cos(qJ(3));
	t69 = t65 * t68;
	t67 = sin(qJ(3));
	t62 = qJ(3) + pkin(9);
	t61 = cos(t62);
	t60 = sin(t62);
	t1 = [0, 0, ((-t60 * t64 - t61 * t70) * r_i_i_C(1) + (t60 * t70 - t61 * t64) * r_i_i_C(2) + (-t64 * t67 - t66 * t69) * pkin(3)) * qJD(3), 0, 0; 0, 0, ((t60 * t66 - t61 * t71) * r_i_i_C(1) + (t60 * t71 + t61 * t66) * r_i_i_C(2) + (-t64 * t69 + t66 * t67) * pkin(3)) * qJD(3), 0, 0; 0, 0, (-pkin(3) * t68 - r_i_i_C(1) * t61 + r_i_i_C(2) * t60) * sin(pkin(8)) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:42
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (115->32), mult. (208->65), div. (0->0), fcn. (184->10), ass. (0->31)
	t190 = qJ(3) + pkin(9);
	t189 = cos(t190);
	t192 = sin(pkin(7));
	t193 = cos(pkin(8));
	t188 = sin(t190);
	t194 = cos(pkin(7));
	t207 = t194 * t188;
	t215 = -t192 * t189 + t193 * t207;
	t214 = -pkin(6) - r_i_i_C(3);
	t213 = pkin(3) * qJD(3);
	t191 = sin(pkin(8));
	t195 = sin(qJ(5));
	t212 = t191 * t195;
	t197 = cos(qJ(5));
	t211 = t191 * t197;
	t209 = t192 * t193;
	t198 = cos(qJ(3));
	t208 = t193 * t198;
	t206 = t194 * t189;
	t204 = t195 * r_i_i_C(1) + t197 * r_i_i_C(2);
	t203 = t197 * r_i_i_C(1) - t195 * r_i_i_C(2) + pkin(4);
	t202 = t204 * t188;
	t186 = t192 * t188 + t193 * t206;
	t184 = t189 * t209 - t207;
	t201 = t188 * t209 + t206;
	t200 = qJD(5) * t204;
	t199 = qJD(3) * t203;
	t196 = sin(qJ(3));
	t181 = t215 * qJD(3);
	t179 = t201 * qJD(3);
	t1 = [0, 0, t214 * t181 + t215 * t200 + (-t192 * t196 - t194 * t208) * t213 - t186 * t199, 0, t204 * t181 + ((-t186 * t197 - t194 * t212) * r_i_i_C(1) + (t186 * t195 - t194 * t211) * r_i_i_C(2)) * qJD(5); 0, 0, t214 * t179 + t201 * t200 + (-t192 * t208 + t194 * t196) * t213 - t184 * t199, 0, t204 * t179 + ((-t184 * t197 - t192 * t212) * r_i_i_C(1) + (t184 * t195 - t192 * t211) * r_i_i_C(2)) * qJD(5); 0, 0, (qJD(5) * t202 + (-pkin(3) * t198 + t214 * t188 - t203 * t189) * qJD(3)) * t191, 0, t191 * qJD(3) * t202 + ((-t189 * t211 + t193 * t195) * r_i_i_C(1) + (t189 * t212 + t193 * t197) * r_i_i_C(2)) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end