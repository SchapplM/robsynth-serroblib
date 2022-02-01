% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:35:07
	% EndTime: 2022-01-23 09:35:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:35:07
	% EndTime: 2022-01-23 09:35:07
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
	% StartTime: 2022-01-23 09:35:07
	% EndTime: 2022-01-23 09:35:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(9);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:35:07
	% EndTime: 2022-01-23 09:35:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (34->9), mult. (24->12), div. (0->0), fcn. (12->6), ass. (0->8)
	t42 = qJ(1) + pkin(9);
	t40 = qJ(3) + t42;
	t38 = sin(t40);
	t39 = cos(t40);
	t41 = qJD(1) + qJD(3);
	t44 = (-r_i_i_C(1) * t39 + r_i_i_C(2) * t38) * t41;
	t43 = (-r_i_i_C(1) * t38 - r_i_i_C(2) * t39) * t41;
	t1 = [(-cos(t42) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t44, 0, t44, 0, 0; t43 + (-sin(t42) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1), 0, t43, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:35:07
	% EndTime: 2022-01-23 09:35:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (88->13), mult. (40->15), div. (0->0), fcn. (20->8), ass. (0->13)
	t50 = qJD(1) + qJD(3);
	t56 = pkin(3) * t50;
	t51 = qJ(1) + pkin(9);
	t49 = qJ(3) + t51;
	t47 = qJ(4) + t49;
	t43 = sin(t47);
	t44 = cos(t47);
	t48 = qJD(4) + t50;
	t55 = (-r_i_i_C(1) * t44 + r_i_i_C(2) * t43) * t48;
	t54 = (-r_i_i_C(1) * t43 - r_i_i_C(2) * t44) * t48;
	t53 = -cos(t49) * t56 + t55;
	t52 = -sin(t49) * t56 + t54;
	t1 = [(-cos(t51) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t53, 0, t53, t55, 0; (-sin(t51) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1) + t52, 0, t52, t54, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:35:07
	% EndTime: 2022-01-23 09:35:08
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (283->26), mult. (152->39), div. (0->0), fcn. (92->10), ass. (0->24)
	t74 = r_i_i_C(3) + pkin(8);
	t53 = qJ(1) + pkin(9);
	t51 = qJ(3) + t53;
	t49 = qJ(4) + t51;
	t45 = sin(t49);
	t46 = cos(t49);
	t55 = cos(qJ(5));
	t66 = qJD(5) * t55;
	t52 = qJD(1) + qJD(3);
	t50 = qJD(4) + t52;
	t54 = sin(qJ(5));
	t69 = t50 * t54;
	t73 = t45 * t69 - t46 * t66;
	t72 = t45 * t66 + t46 * t69;
	t71 = pkin(3) * t52;
	t68 = t50 * t55;
	t67 = qJD(5) * t54;
	t63 = t45 * t67;
	t60 = t45 * t68 + t46 * t67;
	t59 = r_i_i_C(1) * t63 + ((-r_i_i_C(1) * t55 - pkin(4)) * t46 - t74 * t45) * t50 + t72 * r_i_i_C(2);
	t58 = -cos(t51) * t71 + t59;
	t57 = -t60 * r_i_i_C(1) + t73 * r_i_i_C(2) + (-pkin(4) * t45 + t46 * t74) * t50;
	t56 = -sin(t51) * t71 + t57;
	t1 = [(-cos(t53) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t58, 0, t58, t59, t73 * r_i_i_C(1) + t60 * r_i_i_C(2); (-sin(t53) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1) + t56, 0, t56, t57, (-t46 * t68 + t63) * r_i_i_C(2) - t72 * r_i_i_C(1); 0, 0, 0, 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * qJD(5);];
	JaD_transl = t1;
end