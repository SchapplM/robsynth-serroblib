% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:51
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPPR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPPR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:51:26
	% EndTime: 2019-12-29 15:51:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:51:26
	% EndTime: 2019-12-29 15:51:26
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
	% StartTime: 2019-12-29 15:51:26
	% EndTime: 2019-12-29 15:51:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) - r_i_i_C(1);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:51:31
	% EndTime: 2019-12-29 15:51:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (14->11), mult. (36->20), div. (0->0), fcn. (26->4), ass. (0->8)
	t58 = -pkin(1) - pkin(2);
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t55 = cos(pkin(7));
	t54 = sin(pkin(7));
	t53 = (t54 * t57 - t55 * t56) * qJD(1);
	t52 = (t54 * t56 + t55 * t57) * qJD(1);
	t1 = [-t52 * r_i_i_C(1) + t53 * r_i_i_C(2) + t57 * qJD(2) + (-qJ(2) * t56 + t57 * t58) * qJD(1), qJD(1) * t57, 0, 0, 0; t53 * r_i_i_C(1) + t52 * r_i_i_C(2) + t56 * qJD(2) + (qJ(2) * t57 + t56 * t58) * qJD(1), qJD(1) * t56, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:51:26
	% EndTime: 2019-12-29 15:51:26
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (32->17), mult. (84->24), div. (0->0), fcn. (70->6), ass. (0->12)
	t36 = -pkin(1) - pkin(2);
	t35 = -r_i_i_C(3) - qJ(4);
	t27 = sin(pkin(7));
	t29 = cos(pkin(7));
	t30 = sin(qJ(1));
	t31 = cos(qJ(1));
	t34 = t31 * t27 - t30 * t29;
	t33 = t30 * t27 + t31 * t29;
	t32 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(3);
	t25 = t34 * qJD(1);
	t24 = t33 * qJD(1);
	t1 = [-t33 * qJD(4) + t31 * qJD(2) + t35 * t25 - t32 * t24 + (-t30 * qJ(2) + t36 * t31) * qJD(1), qJD(1) * t31, 0, -t24, 0; t34 * qJD(4) + t30 * qJD(2) + t35 * t24 + t32 * t25 + (t31 * qJ(2) + t36 * t30) * qJD(1), qJD(1) * t30, 0, t25, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:51:31
	% EndTime: 2019-12-29 15:51:31
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (77->28), mult. (154->44), div. (0->0), fcn. (134->7), ass. (0->18)
	t52 = pkin(8) + qJ(5);
	t50 = sin(t52);
	t51 = cos(t52);
	t58 = (r_i_i_C(1) * t50 + r_i_i_C(2) * t51) * qJD(5);
	t64 = -pkin(1) - pkin(2);
	t63 = -r_i_i_C(3) - pkin(6) - qJ(4);
	t62 = qJD(5) * t50;
	t61 = qJD(5) * t51;
	t53 = sin(pkin(7));
	t54 = cos(pkin(7));
	t56 = sin(qJ(1));
	t57 = cos(qJ(1));
	t46 = t53 * t57 - t56 * t54;
	t47 = t56 * t53 + t54 * t57;
	t59 = r_i_i_C(1) * t51 - r_i_i_C(2) * t50 + cos(pkin(8)) * pkin(4) + pkin(3);
	t45 = t46 * qJD(1);
	t44 = t47 * qJD(1);
	t1 = [qJD(2) * t57 - t47 * qJD(4) + t63 * t45 - t46 * t58 - t59 * t44 + (-qJ(2) * t56 + t64 * t57) * qJD(1), qJD(1) * t57, 0, -t44, (-t45 * t51 + t47 * t62) * r_i_i_C(2) + (-t45 * t50 - t47 * t61) * r_i_i_C(1); t56 * qJD(2) + t46 * qJD(4) + t63 * t44 - t47 * t58 + t59 * t45 + (qJ(2) * t57 + t64 * t56) * qJD(1), qJD(1) * t56, 0, t45, (-t44 * t51 - t46 * t62) * r_i_i_C(2) + (-t44 * t50 + t46 * t61) * r_i_i_C(1); 0, 0, 0, 0, t58;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end