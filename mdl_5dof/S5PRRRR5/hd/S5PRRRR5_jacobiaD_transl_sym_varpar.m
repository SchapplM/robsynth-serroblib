% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR5
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:36
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR5_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:57
	% EndTime: 2019-10-24 10:36:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:57
	% EndTime: 2019-10-24 10:36:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:57
	% EndTime: 2019-10-24 10:36:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->3), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->4)
	t32 = pkin(9) + qJ(2);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [0, (-r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(2), 0, 0, 0; 0, (-r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(2), 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:57
	% EndTime: 2019-10-24 10:36:58
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (41->16), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t24 = sin(qJ(3));
	t25 = cos(qJ(3));
	t26 = (r_i_i_C(1) * t24 + r_i_i_C(2) * t25) * qJD(3);
	t33 = pkin(6) + r_i_i_C(3);
	t32 = qJD(2) * t24;
	t31 = qJD(2) * t25;
	t30 = qJD(3) * t24;
	t29 = qJD(3) * t25;
	t27 = -r_i_i_C(1) * t25 + r_i_i_C(2) * t24 - pkin(2);
	t23 = pkin(9) + qJ(2);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [0, t21 * t26 + (-t33 * t21 + t27 * t22) * qJD(2), (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0; 0, -t22 * t26 + (t27 * t21 + t33 * t22) * qJD(2), (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0; 0, 0, -t26, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:57
	% EndTime: 2019-10-24 10:36:57
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (117->27), mult. (114->38), div. (0->0), fcn. (73->6), ass. (0->26)
	t45 = qJD(3) + qJD(4);
	t46 = qJ(3) + qJ(4);
	t42 = sin(t46);
	t43 = cos(t46);
	t63 = r_i_i_C(2) * t43;
	t54 = r_i_i_C(1) * t42 + t63;
	t52 = t54 * t45;
	t47 = sin(qJ(3));
	t65 = pkin(3) * t47;
	t66 = qJD(3) * t65 + t52;
	t64 = r_i_i_C(2) * t42;
	t62 = r_i_i_C(3) + pkin(7) + pkin(6);
	t61 = t43 * t45;
	t60 = qJD(2) * t42;
	t48 = cos(qJ(3));
	t59 = qJD(3) * t48;
	t58 = r_i_i_C(1) * t61;
	t57 = t45 * t64;
	t56 = qJD(2) * t63;
	t53 = -t48 * pkin(3) - r_i_i_C(1) * t43 - pkin(2) + t64;
	t44 = pkin(9) + qJ(2);
	t40 = sin(t44);
	t41 = cos(t44);
	t51 = (t57 - t58) * t41 + (t60 * r_i_i_C(1) + t56) * t40;
	t35 = t40 * t57;
	t1 = [0, t66 * t40 + (-t62 * t40 + t53 * t41) * qJD(2), (qJD(2) * t40 * t47 - t41 * t59) * pkin(3) + t51, t51, 0; 0, -t66 * t41 + (t53 * t40 + t62 * t41) * qJD(2), t35 + (-pkin(3) * t59 - t58) * t40 + (-t54 - t65) * t41 * qJD(2), -t41 * t56 + t35 + (-t40 * t61 - t41 * t60) * r_i_i_C(1), 0; 0, 0, -t66, -t52, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:57
	% EndTime: 2019-10-24 10:36:58
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (250->37), mult. (170->47), div. (0->0), fcn. (109->8), ass. (0->38)
	t61 = sin(qJ(3));
	t58 = qJD(3) + qJD(4);
	t53 = qJD(5) + t58;
	t60 = qJ(3) + qJ(4);
	t56 = qJ(5) + t60;
	t50 = cos(t56);
	t82 = r_i_i_C(2) * t50;
	t49 = sin(t56);
	t84 = r_i_i_C(1) * t49;
	t67 = t82 + t84;
	t65 = t67 * t53;
	t54 = sin(t60);
	t85 = pkin(4) * t54;
	t63 = -t58 * t85 - t65;
	t78 = pkin(3) * qJD(3);
	t87 = -t61 * t78 + t63;
	t80 = t50 * t53;
	t73 = r_i_i_C(1) * t80;
	t55 = cos(t60);
	t79 = t55 * t58;
	t86 = -pkin(4) * t79 - t73;
	t83 = r_i_i_C(2) * t49;
	t81 = r_i_i_C(3) + pkin(8) + pkin(7) + pkin(6);
	t57 = pkin(9) + qJ(2);
	t51 = sin(t57);
	t77 = qJD(2) * t51;
	t52 = cos(t57);
	t76 = qJD(2) * t52;
	t72 = t53 * t83;
	t70 = qJD(2) * t82;
	t71 = t51 * t70 + t52 * t72 + t77 * t84;
	t62 = cos(qJ(3));
	t68 = -t62 * t78 + t86;
	t66 = -pkin(3) * t62 - pkin(4) * t55 - r_i_i_C(1) * t50 - pkin(2) + t83;
	t64 = -t52 * t73 + t71;
	t48 = -pkin(3) * t61 - t85;
	t41 = t51 * t72;
	t1 = [0, -t87 * t51 + (-t81 * t51 + t66 * t52) * qJD(2), -t48 * t77 + t68 * t52 + t71, (-t52 * t79 + t54 * t77) * pkin(4) + t64, t64; 0, t87 * t52 + (t66 * t51 + t81 * t52) * qJD(2), t41 + t68 * t51 + (t48 - t67) * t76, t41 + t86 * t51 + (-t67 - t85) * t76, -t52 * t70 + t41 + (-t49 * t76 - t51 * t80) * r_i_i_C(1); 0, 0, t87, t63, -t65;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end