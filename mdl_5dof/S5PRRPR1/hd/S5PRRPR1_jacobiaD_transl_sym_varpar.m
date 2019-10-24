% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPR1
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
% Datum: 2019-10-24 10:29
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:54
	% EndTime: 2019-10-24 10:29:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:54
	% EndTime: 2019-10-24 10:29:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:54
	% EndTime: 2019-10-24 10:29:54
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
	% StartTime: 2019-10-24 10:29:54
	% EndTime: 2019-10-24 10:29:54
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (32->7), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->9)
	t45 = pkin(2) * qJD(2);
	t41 = pkin(8) + qJ(2);
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
	% StartTime: 2019-10-24 10:29:54
	% EndTime: 2019-10-24 10:29:54
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (98->12), mult. (58->14), div. (0->0), fcn. (36->6), ass. (0->14)
	t38 = pkin(8) + qJ(2);
	t37 = qJ(3) + t38;
	t35 = sin(t37);
	t39 = qJD(2) + qJD(3);
	t48 = t39 * t35;
	t50 = qJ(4) + r_i_i_C(3);
	t49 = qJD(4) + r_i_i_C(2) * t39 * sin(pkin(9));
	t36 = cos(t37);
	t47 = t39 * t36;
	t46 = pkin(2) * qJD(2);
	t44 = -r_i_i_C(1) * cos(pkin(9)) - pkin(3);
	t43 = t35 * t49 + t44 * t48 + t47 * t50;
	t42 = -t50 * t48 + (t39 * t44 + t49) * t36;
	t1 = [0, -cos(t38) * t46 + t42, t42, t47, 0; 0, -sin(t38) * t46 + t43, t43, t48, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:54
	% EndTime: 2019-10-24 10:29:54
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (183->26), mult. (114->37), div. (0->0), fcn. (74->7), ass. (0->22)
	t51 = pkin(8) + qJ(2);
	t49 = qJ(3) + t51;
	t44 = sin(t49);
	t50 = pkin(9) + qJ(5);
	t47 = sin(t50);
	t48 = cos(t50);
	t62 = qJD(5) * t48;
	t45 = cos(t49);
	t52 = qJD(2) + qJD(3);
	t66 = t52 * t45;
	t68 = t44 * t62 + t47 * t66;
	t67 = t52 * t44;
	t65 = t52 * (-pkin(7) - qJ(4));
	t64 = pkin(2) * qJD(2);
	t63 = qJD(5) * t47;
	t61 = t47 * t67;
	t59 = t44 * t63;
	t57 = -r_i_i_C(1) * t48 - cos(pkin(9)) * pkin(4) - pkin(3);
	t56 = (-r_i_i_C(1) * t47 - r_i_i_C(2) * t48) * qJD(5);
	t55 = r_i_i_C(1) * t59 + t44 * t65 + t45 * qJD(4) + (-r_i_i_C(3) * t44 + t45 * t57) * t52 + t68 * r_i_i_C(2);
	t54 = r_i_i_C(2) * t61 + r_i_i_C(3) * t66 + t44 * qJD(4) + t57 * t67 + (t56 - t65) * t45;
	t1 = [0, -cos(t51) * t64 + t55, t55, t66, (t45 * t63 + t48 * t67) * r_i_i_C(2) + (-t45 * t62 + t61) * r_i_i_C(1); 0, -sin(t51) * t64 + t54, t54, t67, (-t48 * t66 + t59) * r_i_i_C(2) - t68 * r_i_i_C(1); 0, 0, 0, 0, t56;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end