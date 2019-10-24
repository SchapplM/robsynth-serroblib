% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR3
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

function JaD_transl = S5PRRRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:13
	% EndTime: 2019-10-24 10:36:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:13
	% EndTime: 2019-10-24 10:36:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:13
	% EndTime: 2019-10-24 10:36:13
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
	% StartTime: 2019-10-24 10:36:13
	% EndTime: 2019-10-24 10:36:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (32->7), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->9)
	t45 = pkin(2) * qJD(2);
	t41 = pkin(9) + qJ(2);
	t40 = qJ(3) + t41;
	t38 = sin(t40);
	t39 = cos(t40);
	t42 = qJD(2) + qJD(3);
	t44 = (-r_i_i_C(1) * t39 + r_i_i_C(2) * t38) * t42;
	t43 = (-r_i_i_C(1) * t38 - r_i_i_C(2) * t39) * t42;
	t1 = [0, -cos(t41) * t45 + t44, t44, 0, 0; 0, -sin(t41) * t45 + t43, t43, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:13
	% EndTime: 2019-10-24 10:36:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (86->11), mult. (36->12), div. (0->0), fcn. (18->6), ass. (0->14)
	t51 = qJD(2) + qJD(3);
	t57 = pkin(3) * t51;
	t56 = pkin(2) * qJD(2);
	t50 = pkin(9) + qJ(2);
	t49 = qJ(3) + t50;
	t47 = qJ(4) + t49;
	t43 = sin(t47);
	t44 = cos(t47);
	t48 = qJD(4) + t51;
	t55 = (-r_i_i_C(1) * t44 + r_i_i_C(2) * t43) * t48;
	t54 = (-r_i_i_C(1) * t43 - r_i_i_C(2) * t44) * t48;
	t53 = -cos(t49) * t57 + t55;
	t52 = -sin(t49) * t57 + t54;
	t1 = [0, -cos(t50) * t56 + t53, t53, t55, 0; 0, -sin(t50) * t56 + t52, t52, t54, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:13
	% EndTime: 2019-10-24 10:36:13
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (281->24), mult. (148->36), div. (0->0), fcn. (90->8), ass. (0->25)
	t75 = r_i_i_C(3) + pkin(8);
	t52 = pkin(9) + qJ(2);
	t51 = qJ(3) + t52;
	t49 = qJ(4) + t51;
	t45 = sin(t49);
	t46 = cos(t49);
	t55 = cos(qJ(5));
	t66 = qJD(5) * t55;
	t53 = qJD(2) + qJD(3);
	t50 = qJD(4) + t53;
	t54 = sin(qJ(5));
	t70 = t50 * t54;
	t74 = t45 * t70 - t46 * t66;
	t73 = t45 * t66 + t46 * t70;
	t72 = pkin(3) * t53;
	t69 = t50 * t55;
	t68 = pkin(2) * qJD(2);
	t67 = qJD(5) * t54;
	t63 = t45 * t67;
	t60 = t45 * t69 + t46 * t67;
	t59 = r_i_i_C(1) * t63 + ((-r_i_i_C(1) * t55 - pkin(4)) * t46 - t75 * t45) * t50 + t73 * r_i_i_C(2);
	t58 = -cos(t51) * t72 + t59;
	t57 = -t60 * r_i_i_C(1) + t74 * r_i_i_C(2) + (-pkin(4) * t45 + t46 * t75) * t50;
	t56 = -sin(t51) * t72 + t57;
	t1 = [0, -cos(t52) * t68 + t58, t58, t59, t74 * r_i_i_C(1) + t60 * r_i_i_C(2); 0, -sin(t52) * t68 + t56, t56, t57, (-t46 * t69 + t63) * r_i_i_C(2) - t73 * r_i_i_C(1); 0, 0, 0, 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end