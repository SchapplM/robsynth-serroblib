% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRP3
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:44:39
	% EndTime: 2019-12-05 16:44:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:44:39
	% EndTime: 2019-12-05 16:44:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:44:39
	% EndTime: 2019-12-05 16:44:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->3), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->4)
	t32 = pkin(8) + qJ(2);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [0, (-r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(2), 0, 0, 0; 0, (-r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(2), 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:44:39
	% EndTime: 2019-12-05 16:44:39
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
	t23 = pkin(8) + qJ(2);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [0, t21 * t26 + (-t21 * t33 + t22 * t27) * qJD(2), (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0; 0, -t22 * t26 + (t21 * t27 + t22 * t33) * qJD(2), (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0; 0, 0, -t26, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:44:39
	% EndTime: 2019-12-05 16:44:39
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
	t44 = pkin(8) + qJ(2);
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
	% StartTime: 2019-12-05 16:44:39
	% EndTime: 2019-12-05 16:44:39
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (165->33), mult. (146->44), div. (0->0), fcn. (95->6), ass. (0->32)
	t53 = qJ(3) + qJ(4);
	t49 = cos(t53);
	t52 = qJD(3) + qJD(4);
	t69 = t49 * t52;
	t73 = -pkin(4) - r_i_i_C(1);
	t75 = t73 * t69;
	t48 = sin(t53);
	t72 = r_i_i_C(2) * t49;
	t59 = r_i_i_C(1) * t48 + t72;
	t54 = sin(qJ(3));
	t68 = pkin(3) * qJD(3);
	t62 = t54 * t68;
	t70 = t48 * t52;
	t74 = pkin(4) * t70 + t59 * t52 + t62;
	t71 = r_i_i_C(3) + qJ(5) + pkin(7) + pkin(6);
	t51 = pkin(8) + qJ(2);
	t46 = sin(t51);
	t67 = qJD(2) * t46;
	t47 = cos(t51);
	t66 = qJD(2) * t47;
	t65 = r_i_i_C(2) * t70;
	t64 = t47 * t69;
	t61 = t48 * t67;
	t63 = r_i_i_C(1) * t61 + t47 * t65 + t67 * t72;
	t55 = cos(qJ(3));
	t60 = -t55 * t68 + t75;
	t58 = -pkin(3) * t55 + r_i_i_C(2) * t48 + t73 * t49 - pkin(2);
	t57 = t73 * t48 - t72;
	t56 = t57 * t52;
	t45 = -pkin(3) * t54 - pkin(4) * t48;
	t40 = t46 * t65;
	t1 = [0, t47 * qJD(5) + t74 * t46 + (-t71 * t46 + t58 * t47) * qJD(2), -t45 * t67 + t60 * t47 + t63, -r_i_i_C(1) * t64 + (t61 - t64) * pkin(4) + t63, t66; 0, t46 * qJD(5) - t74 * t47 + (t58 * t46 + t71 * t47) * qJD(2), t40 + t60 * t46 + (t45 - t59) * t66, t46 * t75 + t57 * t66 + t40, t67; 0, 0, t56 - t62, t56, 0;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end