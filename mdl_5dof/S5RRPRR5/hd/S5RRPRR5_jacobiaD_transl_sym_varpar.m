% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:03:43
	% EndTime: 2022-01-20 11:03:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:03:43
	% EndTime: 2022-01-20 11:03:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:03:43
	% EndTime: 2022-01-20 11:03:43
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t40 = qJ(1) + qJ(2);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(1) + qJD(2);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [-cos(qJ(1)) * t43 + t42, t42, 0, 0, 0; -sin(qJ(1)) * t43 + t41, t41, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:03:43
	% EndTime: 2022-01-20 11:03:43
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (70->11), mult. (58->14), div. (0->0), fcn. (36->6), ass. (0->13)
	t37 = qJ(1) + qJ(2);
	t34 = sin(t37);
	t36 = qJD(1) + qJD(2);
	t46 = t36 * t34;
	t48 = qJ(3) + r_i_i_C(3);
	t47 = qJD(3) + r_i_i_C(2) * t36 * sin(pkin(9));
	t35 = cos(t37);
	t45 = t36 * t35;
	t44 = pkin(1) * qJD(1);
	t42 = -r_i_i_C(1) * cos(pkin(9)) - pkin(2);
	t41 = t47 * t34 + t42 * t46 + t48 * t45;
	t40 = -t48 * t46 + (t42 * t36 + t47) * t35;
	t1 = [-cos(qJ(1)) * t44 + t40, t40, t45, 0, 0; -sin(qJ(1)) * t44 + t41, t41, t46, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:03:43
	% EndTime: 2022-01-20 11:03:43
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (139->25), mult. (114->37), div. (0->0), fcn. (74->7), ass. (0->21)
	t48 = pkin(9) + qJ(4);
	t44 = sin(t48);
	t45 = cos(t48);
	t50 = qJ(1) + qJ(2);
	t46 = sin(t50);
	t61 = qJD(4) * t46;
	t47 = cos(t50);
	t49 = qJD(1) + qJD(2);
	t64 = t49 * t47;
	t66 = t44 * t64 + t45 * t61;
	t65 = t49 * t46;
	t63 = t49 * (-pkin(7) - qJ(3));
	t62 = pkin(1) * qJD(1);
	t60 = qJD(4) * t47;
	t59 = t44 * t65;
	t57 = t44 * t61;
	t55 = -r_i_i_C(1) * t45 - cos(pkin(9)) * pkin(3) - pkin(2);
	t54 = (-r_i_i_C(1) * t44 - r_i_i_C(2) * t45) * qJD(4);
	t53 = r_i_i_C(1) * t57 + t46 * t63 + t47 * qJD(3) + (-r_i_i_C(3) * t46 + t55 * t47) * t49 + t66 * r_i_i_C(2);
	t52 = r_i_i_C(2) * t59 + r_i_i_C(3) * t64 + t46 * qJD(3) + t55 * t65 + (t54 - t63) * t47;
	t1 = [-cos(qJ(1)) * t62 + t53, t53, t64, (t44 * t60 + t45 * t65) * r_i_i_C(2) + (-t45 * t60 + t59) * r_i_i_C(1), 0; -sin(qJ(1)) * t62 + t52, t52, t65, (-t45 * t64 + t57) * r_i_i_C(2) - t66 * r_i_i_C(1), 0; 0, 0, 0, t54, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:03:43
	% EndTime: 2022-01-20 11:03:43
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (264->36), mult. (172->48), div. (0->0), fcn. (113->9), ass. (0->35)
	t66 = pkin(9) + qJ(4);
	t62 = qJ(5) + t66;
	t58 = sin(t62);
	t59 = cos(t62);
	t69 = qJ(1) + qJ(2);
	t64 = cos(t69);
	t68 = qJD(1) + qJD(2);
	t86 = t68 * t64;
	t63 = sin(t69);
	t67 = qJD(4) + qJD(5);
	t88 = t63 * t67;
	t92 = t58 * t86 + t59 * t88;
	t60 = sin(t66);
	t91 = pkin(4) * t60;
	t90 = r_i_i_C(1) * t59;
	t89 = r_i_i_C(2) * t59;
	t87 = t68 * t63;
	t85 = pkin(1) * qJD(1);
	t61 = cos(t66);
	t84 = qJD(4) * t61;
	t83 = t67 * t90;
	t82 = t68 * t89;
	t81 = t58 * t88;
	t80 = t58 * t87;
	t77 = qJD(4) * t91;
	t76 = -pkin(4) * t61 - cos(pkin(9)) * pkin(3) - pkin(2) - t90;
	t75 = -r_i_i_C(1) * t58 - t89;
	t74 = t75 * t67;
	t73 = r_i_i_C(1) * t80 + t63 * t82 + (r_i_i_C(2) * t58 * t67 - t83) * t64;
	t72 = t74 - t77;
	t65 = qJ(3) + pkin(7) + pkin(8);
	t71 = r_i_i_C(1) * t81 + t63 * t77 + qJD(3) * t64 + (t76 * t64 + (-r_i_i_C(3) - t65) * t63) * t68 + t92 * r_i_i_C(2);
	t70 = r_i_i_C(2) * t80 + r_i_i_C(3) * t86 + t63 * qJD(3) + t76 * t87 + (t65 * t68 + t72) * t64;
	t45 = r_i_i_C(2) * t81;
	t1 = [-cos(qJ(1)) * t85 + t71, t71, t86, (t60 * t87 - t64 * t84) * pkin(4) + t73, t73; -sin(qJ(1)) * t85 + t70, t70, t87, t45 + (-pkin(4) * t84 - t83) * t63 + (t75 - t91) * t86, -r_i_i_C(1) * t92 - t64 * t82 + t45; 0, 0, 0, t72, t74;];
	JaD_transl = t1;
end