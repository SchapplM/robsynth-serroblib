% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:30
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:38
	% EndTime: 2019-10-24 10:30:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:38
	% EndTime: 2019-10-24 10:30:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:38
	% EndTime: 2019-10-24 10:30:38
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
	% StartTime: 2019-10-24 10:30:38
	% EndTime: 2019-10-24 10:30:38
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
	t1 = [0, t21 * t26 + (-t33 * t21 + t27 * t22) * qJD(2), (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0; 0, -t22 * t26 + (t27 * t21 + t33 * t22) * qJD(2), (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0; 0, 0, -t26, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:38
	% EndTime: 2019-10-24 10:30:38
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (82->17), mult. (90->25), div. (0->0), fcn. (59->6), ass. (0->15)
	t34 = qJ(3) + pkin(9);
	t30 = sin(t34);
	t32 = cos(t34);
	t47 = -r_i_i_C(1) * t32 + r_i_i_C(2) * t30 - cos(qJ(3)) * pkin(3);
	t45 = r_i_i_C(3) + qJ(4) + pkin(6);
	t33 = pkin(8) + qJ(2);
	t31 = cos(t33);
	t44 = qJD(2) * t31;
	t42 = -pkin(2) + t47;
	t41 = sin(qJ(3)) * pkin(3) + r_i_i_C(1) * t30 + r_i_i_C(2) * t32;
	t29 = sin(t33);
	t40 = t41 * t29;
	t39 = qJD(3) * t47;
	t38 = t41 * qJD(3);
	t1 = [0, t31 * qJD(4) + qJD(3) * t40 + (-t29 * t45 + t31 * t42) * qJD(2), qJD(2) * t40 + t31 * t39, t44, 0; 0, t29 * qJD(4) - t31 * t38 + (t29 * t42 + t31 * t45) * qJD(2), t29 * t39 - t41 * t44, qJD(2) * t29, 0; 0, 0, -t38, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:38
	% EndTime: 2019-10-24 10:30:38
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (171->31), mult. (132->40), div. (0->0), fcn. (86->8), ass. (0->28)
	t56 = qJ(3) + pkin(9);
	t45 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t56);
	t55 = qJD(3) + qJD(5);
	t52 = qJ(5) + t56;
	t47 = cos(t52);
	t72 = r_i_i_C(2) * t47;
	t46 = sin(t52);
	t74 = r_i_i_C(1) * t46;
	t61 = t72 + t74;
	t59 = t61 * t55;
	t75 = t45 * qJD(3) - t59;
	t73 = r_i_i_C(2) * t46;
	t71 = r_i_i_C(3) + pkin(7) + qJ(4) + pkin(6);
	t70 = t47 * t55;
	t54 = pkin(8) + qJ(2);
	t48 = sin(t54);
	t69 = qJD(2) * t48;
	t50 = cos(t54);
	t68 = qJD(2) * t50;
	t67 = r_i_i_C(1) * t70;
	t66 = t55 * t73;
	t64 = qJD(2) * t72;
	t65 = t48 * t64 + t50 * t66 + t69 * t74;
	t62 = -cos(qJ(3)) * pkin(3) - pkin(4) * cos(t56);
	t63 = t62 * qJD(3) - t67;
	t60 = -r_i_i_C(1) * t47 - pkin(2) + t62 + t73;
	t38 = t48 * t66;
	t1 = [0, t50 * qJD(4) - t75 * t48 + (-t71 * t48 + t60 * t50) * qJD(2), -t45 * t69 + t63 * t50 + t65, t68, -t50 * t67 + t65; 0, t48 * qJD(4) + t75 * t50 + (t60 * t48 + t71 * t50) * qJD(2), t38 + t63 * t48 + (t45 - t61) * t68, t69, -t50 * t64 + t38 + (-t46 * t68 - t48 * t70) * r_i_i_C(1); 0, 0, t75, 0, -t59;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end