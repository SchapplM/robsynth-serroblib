% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:04:53
	% EndTime: 2019-12-29 19:04:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:04:51
	% EndTime: 2019-12-29 19:04:51
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
	% StartTime: 2019-12-29 19:04:53
	% EndTime: 2019-12-29 19:04:53
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-12-29 19:04:58
	% EndTime: 2019-12-29 19:04:58
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (58->11), mult. (42->13), div. (0->0), fcn. (24->4), ass. (0->11)
	t32 = r_i_i_C(3) + qJ(3);
	t26 = qJ(1) + qJ(2);
	t23 = sin(t26);
	t25 = qJD(1) + qJD(2);
	t31 = t25 * t23;
	t24 = cos(t26);
	t30 = t25 * t24;
	t29 = pkin(1) * qJD(1);
	t28 = t23 * qJD(3) + (-pkin(2) + r_i_i_C(2)) * t31 + t32 * t30;
	t27 = r_i_i_C(2) * t30 + t24 * qJD(3) + (-pkin(2) * t24 - t32 * t23) * t25;
	t1 = [-cos(qJ(1)) * t29 + t27, t27, t30, 0, 0; -sin(qJ(1)) * t29 + t28, t28, t31, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:04:51
	% EndTime: 2019-12-29 19:04:52
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (117->23), mult. (118->38), div. (0->0), fcn. (74->6), ass. (0->20)
	t45 = qJ(1) + qJ(2);
	t42 = sin(t45);
	t43 = cos(t45);
	t47 = cos(qJ(4));
	t56 = qJD(4) * t47;
	t44 = qJD(1) + qJD(2);
	t46 = sin(qJ(4));
	t60 = t44 * t46;
	t62 = t42 * t56 + t43 * t60;
	t61 = t44 * t43;
	t59 = t44 * t47;
	t58 = pkin(1) * qJD(1);
	t57 = qJD(4) * t46;
	t55 = -pkin(2) - pkin(7) - r_i_i_C(3);
	t53 = t43 * t59;
	t51 = t43 * t57;
	t50 = t43 * t56;
	t49 = r_i_i_C(2) * t53 + qJ(3) * t61 + (-r_i_i_C(2) * t57 + t55 * t44 + qJD(3)) * t42 + t62 * r_i_i_C(1);
	t48 = -r_i_i_C(2) * t51 + r_i_i_C(1) * t50 + t43 * qJD(3) + (t55 * t43 + (-r_i_i_C(1) * t46 - r_i_i_C(2) * t47 - qJ(3)) * t42) * t44;
	t1 = [-cos(qJ(1)) * t58 + t48, t48, t61, -t62 * r_i_i_C(2) + (-t42 * t57 + t53) * r_i_i_C(1), 0; -sin(qJ(1)) * t58 + t49, t49, t44 * t42, (-t42 * t60 + t50) * r_i_i_C(2) + (t42 * t59 + t51) * r_i_i_C(1), 0; 0, 0, 0, (-r_i_i_C(1) * t47 + r_i_i_C(2) * t46) * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:04:58
	% EndTime: 2019-12-29 19:04:59
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (227->35), mult. (184->46), div. (0->0), fcn. (117->8), ass. (0->34)
	t76 = qJ(1) + qJ(2);
	t70 = sin(t76);
	t74 = qJD(1) + qJD(2);
	t95 = t74 * t70;
	t75 = qJ(4) + qJ(5);
	t69 = sin(t75);
	t72 = cos(t76);
	t94 = t74 * t72;
	t71 = cos(t75);
	t73 = qJD(4) + qJD(5);
	t97 = t71 * t73;
	t105 = t69 * t94 + t70 * t97;
	t77 = sin(qJ(4));
	t104 = pkin(4) * t77 + qJ(3);
	t78 = cos(qJ(4));
	t100 = pkin(4) * t78;
	t85 = qJD(4) * t100;
	t103 = qJD(3) + t85 + (-pkin(8) - pkin(7) - pkin(2) - r_i_i_C(3)) * t74;
	t99 = r_i_i_C(2) * t69;
	t98 = t69 * t73;
	t96 = t72 * t73;
	t92 = pkin(1) * qJD(1);
	t91 = qJD(4) * t77;
	t66 = r_i_i_C(2) * t98;
	t90 = t69 * t96;
	t87 = t71 * t96;
	t86 = t71 * t94;
	t84 = -r_i_i_C(1) * t97 + t66;
	t83 = -r_i_i_C(1) * t69 - r_i_i_C(2) * t71;
	t82 = r_i_i_C(2) * t87 - t95 * t99 + (t71 * t95 + t90) * r_i_i_C(1);
	t81 = r_i_i_C(2) * t86 + t104 * t94 + t105 * r_i_i_C(1) + (t103 - t66) * t70;
	t80 = -r_i_i_C(2) * t90 + r_i_i_C(1) * t87 + (t83 - t104) * t95 + t103 * t72;
	t59 = r_i_i_C(1) * t86;
	t1 = [-cos(qJ(1)) * t92 + t80, t80, t94, t59 + (-t99 + t100) * t94 + (-pkin(4) * t91 + t83 * t73) * t70, -t70 * r_i_i_C(1) * t98 - t105 * r_i_i_C(2) + t59; -sin(qJ(1)) * t92 + t81, t81, t95, (t72 * t91 + t78 * t95) * pkin(4) + t82, t82; 0, 0, 0, t84 - t85, t84;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end