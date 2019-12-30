% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPP4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPP4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:42:50
	% EndTime: 2019-12-29 19:42:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:42:55
	% EndTime: 2019-12-29 19:42:55
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:42:55
	% EndTime: 2019-12-29 19:42:55
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(6) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t18 * t28 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:42:55
	% EndTime: 2019-12-29 19:42:55
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (81->26), mult. (114->37), div. (0->0), fcn. (73->6), ass. (0->27)
	t38 = qJD(2) + qJD(3);
	t39 = qJ(2) + qJ(3);
	t37 = cos(t39);
	t59 = r_i_i_C(2) * t37;
	t36 = sin(t39);
	t61 = r_i_i_C(1) * t36;
	t49 = t59 + t61;
	t47 = t49 * t38;
	t40 = sin(qJ(2));
	t62 = pkin(2) * t40;
	t63 = qJD(2) * t62 + t47;
	t60 = r_i_i_C(2) * t36;
	t58 = r_i_i_C(3) + pkin(7) + pkin(6);
	t57 = t37 * t38;
	t41 = sin(qJ(1));
	t56 = qJD(1) * t41;
	t43 = cos(qJ(1));
	t55 = qJD(1) * t43;
	t42 = cos(qJ(2));
	t54 = qJD(2) * t42;
	t53 = r_i_i_C(1) * t57;
	t52 = t38 * t60;
	t51 = qJD(1) * t59;
	t48 = -t42 * pkin(2) - r_i_i_C(1) * t37 - pkin(1) + t60;
	t46 = t41 * t51 + t56 * t61 + (t52 - t53) * t43;
	t31 = t41 * t52;
	t1 = [t63 * t41 + (-t58 * t41 + t48 * t43) * qJD(1), (t40 * t56 - t43 * t54) * pkin(2) + t46, t46, 0, 0; -t63 * t43 + (t48 * t41 + t58 * t43) * qJD(1), t31 + (-pkin(2) * t54 - t53) * t41 + (-t49 - t62) * t55, -t43 * t51 + t31 + (-t36 * t55 - t41 * t57) * r_i_i_C(1), 0, 0; 0, -t63, -t47, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:42:55
	% EndTime: 2019-12-29 19:42:55
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (149->33), mult. (146->42), div. (0->0), fcn. (95->8), ass. (0->32)
	t52 = sin(qJ(2));
	t50 = qJD(2) + qJD(3);
	t51 = qJ(2) + qJ(3);
	t46 = pkin(8) + t51;
	t44 = sin(t46);
	t45 = cos(t46);
	t59 = r_i_i_C(1) * t44 + r_i_i_C(2) * t45;
	t47 = sin(t51);
	t75 = pkin(3) * t47;
	t78 = t59 + t75;
	t56 = t78 * t50;
	t67 = pkin(2) * qJD(2);
	t76 = t52 * t67 + t56;
	t48 = cos(t51);
	t72 = r_i_i_C(1) * t45;
	t77 = -pkin(3) * t48 - t72;
	t71 = r_i_i_C(2) * t44;
	t69 = r_i_i_C(3) + qJ(4) + pkin(7) + pkin(6);
	t68 = t48 * t50;
	t53 = sin(qJ(1));
	t66 = qJD(1) * t53;
	t55 = cos(qJ(1));
	t65 = qJD(1) * t55;
	t64 = t50 * t72;
	t63 = t50 * t71;
	t62 = t55 * t63 + t59 * t66;
	t54 = cos(qJ(2));
	t60 = -pkin(3) * t68 - t54 * t67 - t64;
	t58 = -t54 * pkin(2) - pkin(1) + t71 + t77;
	t43 = -t52 * pkin(2) - t75;
	t38 = t53 * t63;
	t1 = [t55 * qJD(4) + t76 * t53 + (-t69 * t53 + t58 * t55) * qJD(1), -t43 * t66 + t60 * t55 + t62, -t55 * t64 + (t47 * t66 - t55 * t68) * pkin(3) + t62, t65, 0; t53 * qJD(4) - t76 * t55 + (t58 * t53 + t69 * t55) * qJD(1), t38 + t60 * t53 + (t43 - t59) * t65, t50 * t53 * t77 - t65 * t78 + t38, t66, 0; 0, -t76, -t56, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:42:48
	% EndTime: 2019-12-29 19:42:48
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (288->42), mult. (249->47), div. (0->0), fcn. (170->8), ass. (0->37)
	t195 = qJ(2) + qJ(3);
	t190 = pkin(8) + t195;
	t189 = cos(t190);
	t219 = r_i_i_C(3) + qJ(5);
	t207 = t219 * t189;
	t188 = sin(t190);
	t187 = t188 * qJD(5);
	t191 = sin(t195);
	t194 = qJD(2) + qJD(3);
	t224 = pkin(4) + r_i_i_C(1);
	t212 = t224 * t188;
	t196 = sin(qJ(2));
	t218 = pkin(2) * qJD(2);
	t213 = t196 * t218;
	t221 = pkin(3) * t194;
	t228 = (-t212 + t207) * t194 + (r_i_i_C(2) + qJ(4) + pkin(7) + pkin(6)) * qJD(1) - t191 * t221 + t187 - t213;
	t192 = cos(t195);
	t204 = -t219 * t188 - t224 * t189;
	t226 = -pkin(3) * t192 + t204;
	t223 = pkin(3) * t191;
	t217 = t189 * t194;
	t197 = sin(qJ(1));
	t216 = qJD(1) * t197;
	t199 = cos(qJ(1));
	t215 = qJD(1) * t199;
	t214 = qJD(5) * t189;
	t211 = t197 * t214 + t215 * t207;
	t208 = t188 * t216;
	t210 = t199 * t214 + t224 * t208;
	t205 = -t212 - t223;
	t198 = cos(qJ(2));
	t203 = -t192 * t221 + t204 * t194 - t198 * t218;
	t202 = t205 * t194 + t219 * t217 + t187;
	t201 = t194 * t226;
	t200 = qJD(4) + (-t198 * pkin(2) - pkin(1) + t226) * qJD(1);
	t184 = -t196 * pkin(2) - t223;
	t1 = [-t228 * t197 + t200 * t199, (-t184 - t207) * t216 + t203 * t199 + t210, (-t207 + t223) * t216 + t199 * t201 + t210, t215, t199 * t217 - t208; t200 * t197 + t228 * t199, (t184 - t212) * t215 + t203 * t197 + t211, t197 * t201 + t205 * t215 + t211, t216, t188 * t215 + t197 * t217; 0, t202 - t213, t202, 0, t194 * t188;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end