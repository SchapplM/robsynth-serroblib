% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:52:15
	% EndTime: 2022-01-20 09:52:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:52:15
	% EndTime: 2022-01-20 09:52:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:52:15
	% EndTime: 2022-01-20 09:52:15
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
	% StartTime: 2022-01-20 09:52:14
	% EndTime: 2022-01-20 09:52:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (42->9), mult. (28->11), div. (0->0), fcn. (14->6), ass. (0->9)
	t49 = pkin(1) * qJD(1);
	t46 = qJ(1) + qJ(2);
	t42 = pkin(8) + t46;
	t40 = sin(t42);
	t41 = cos(t42);
	t45 = qJD(1) + qJD(2);
	t48 = (-pkin(2) * cos(t46) - r_i_i_C(1) * t41 + t40 * r_i_i_C(2)) * t45;
	t47 = (-pkin(2) * sin(t46) - r_i_i_C(1) * t40 - r_i_i_C(2) * t41) * t45;
	t1 = [-cos(qJ(1)) * t49 + t48, t48, 0, 0, 0; -sin(qJ(1)) * t49 + t47, t47, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:52:15
	% EndTime: 2022-01-20 09:52:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (108->13), mult. (66->16), div. (0->0), fcn. (40->8), ass. (0->12)
	t53 = qJ(4) + r_i_i_C(3);
	t42 = qJD(1) + qJD(2);
	t52 = qJD(4) + (r_i_i_C(2) * sin(pkin(9)) - pkin(3) - r_i_i_C(1) * cos(pkin(9))) * t42;
	t43 = qJ(1) + qJ(2);
	t39 = pkin(8) + t43;
	t38 = cos(t39);
	t51 = t42 * t38;
	t50 = pkin(1) * qJD(1);
	t37 = sin(t39);
	t47 = -pkin(2) * sin(t43) * t42 + t53 * t51 + t52 * t37;
	t46 = (-pkin(2) * cos(t43) - t53 * t37) * t42 + t52 * t38;
	t1 = [-cos(qJ(1)) * t50 + t46, t46, 0, t51, 0; -sin(qJ(1)) * t50 + t47, t47, 0, t42 * t37, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:52:15
	% EndTime: 2022-01-20 09:52:15
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (193->28), mult. (122->40), div. (0->0), fcn. (78->9), ass. (0->22)
	t56 = qJ(1) + qJ(2);
	t51 = pkin(8) + t56;
	t46 = sin(t51);
	t54 = pkin(9) + qJ(5);
	t49 = sin(t54);
	t50 = cos(t54);
	t66 = qJD(5) * t50;
	t47 = cos(t51);
	t55 = qJD(1) + qJD(2);
	t69 = t55 * t47;
	t71 = t46 * t66 + t49 * t69;
	t70 = t55 * t46;
	t68 = pkin(1) * qJD(1);
	t67 = qJD(5) * t49;
	t65 = t49 * t70;
	t63 = t46 * t67;
	t61 = -r_i_i_C(1) * t50 - cos(pkin(9)) * pkin(4) - pkin(3);
	t60 = (-r_i_i_C(1) * t49 - r_i_i_C(2) * t50) * qJD(5);
	t57 = -pkin(7) - qJ(4);
	t59 = r_i_i_C(1) * t63 + t57 * t70 + t47 * qJD(4) + (-pkin(2) * cos(t56) - r_i_i_C(3) * t46 + t61 * t47) * t55 + t71 * r_i_i_C(2);
	t58 = r_i_i_C(2) * t65 + r_i_i_C(3) * t69 + t46 * qJD(4) + t47 * t60 + (-pkin(2) * sin(t56) + t61 * t46 - t47 * t57) * t55;
	t1 = [-cos(qJ(1)) * t68 + t59, t59, 0, t69, (t47 * t67 + t50 * t70) * r_i_i_C(2) + (-t47 * t66 + t65) * r_i_i_C(1); -sin(qJ(1)) * t68 + t58, t58, 0, t70, (-t50 * t69 + t63) * r_i_i_C(2) - t71 * r_i_i_C(1); 0, 0, 0, 0, t60;];
	JaD_transl = t1;
end