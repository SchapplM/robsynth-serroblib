% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:21:13
	% EndTime: 2020-01-03 11:21:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:21:13
	% EndTime: 2020-01-03 11:21:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (-r_i_i_C(1) * t5 - r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:21:13
	% EndTime: 2020-01-03 11:21:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t9 = qJ(1) + pkin(7);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t7 - r_i_i_C(2) * t8) * qJD(1), 0, 0, 0, 0; (cos(qJ(1)) * pkin(1) + r_i_i_C(1) * t8 - r_i_i_C(2) * t7) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:21:13
	% EndTime: 2020-01-03 11:21:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (27->11), mult. (32->14), div. (0->0), fcn. (20->6), ass. (0->6)
	t42 = r_i_i_C(3) + qJ(3);
	t41 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(2);
	t38 = qJ(1) + pkin(7);
	t37 = cos(t38);
	t36 = sin(t38);
	t1 = [0, 0, 0, 0, 0; t36 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t42 * t37 - t41 * t36) * qJD(1), 0, qJD(1) * t36, 0, 0; -t37 * qJD(3) + (cos(qJ(1)) * pkin(1) + t42 * t36 + t41 * t37) * qJD(1), 0, -qJD(1) * t37, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:21:13
	% EndTime: 2020-01-03 11:21:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (49->17), mult. (66->24), div. (0->0), fcn. (50->8), ass. (0->11)
	t76 = sin(pkin(9));
	t77 = sin(pkin(8));
	t78 = cos(pkin(9));
	t85 = (r_i_i_C(3) + qJ(4)) * t77 + (r_i_i_C(1) * t78 - r_i_i_C(2) * t76 + pkin(3)) * cos(pkin(8)) + pkin(2);
	t83 = qJD(1) * t77;
	t82 = t77 * qJD(4);
	t80 = t76 * r_i_i_C(1) + t78 * r_i_i_C(2) + qJ(3);
	t75 = qJ(1) + pkin(7);
	t74 = cos(t75);
	t73 = sin(t75);
	t1 = [0, 0, 0, 0, 0; t74 * t82 + t73 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t80 * t74 - t85 * t73) * qJD(1), 0, qJD(1) * t73, t74 * t83, 0; t73 * t82 - t74 * qJD(3) + (cos(qJ(1)) * pkin(1) + t80 * t73 + t85 * t74) * qJD(1), 0, -qJD(1) * t74, t73 * t83, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:21:13
	% EndTime: 2020-01-03 11:21:13
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (138->32), mult. (144->52), div. (0->0), fcn. (120->10), ass. (0->23)
	t128 = sin(pkin(8));
	t129 = cos(pkin(8));
	t142 = (r_i_i_C(3) + pkin(6) + qJ(4)) * t128 + (cos(pkin(9)) * pkin(4) + pkin(3)) * t129 + pkin(2);
	t126 = qJ(1) + pkin(7);
	t122 = sin(t126);
	t140 = t122 * t129;
	t124 = cos(t126);
	t139 = t124 * t129;
	t138 = qJD(1) * t128;
	t137 = qJD(4) * t128;
	t135 = pkin(4) * sin(pkin(9)) + qJ(3);
	t125 = pkin(9) + qJ(5);
	t121 = sin(t125);
	t123 = cos(t125);
	t134 = t121 * t122 + t123 * t139;
	t133 = t121 * t124 - t123 * t140;
	t132 = -t121 * t139 + t122 * t123;
	t131 = t121 * t140 + t123 * t124;
	t119 = t134 * qJD(1) - t131 * qJD(5);
	t118 = t132 * qJD(1) + t133 * qJD(5);
	t117 = t133 * qJD(1) + t132 * qJD(5);
	t116 = t131 * qJD(1) - t134 * qJD(5);
	t1 = [0, 0, 0, 0, (-r_i_i_C(1) * t123 + r_i_i_C(2) * t121) * t128 * qJD(5); t124 * t137 + t117 * r_i_i_C(1) + t116 * r_i_i_C(2) + t122 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t135 * t124 - t142 * t122) * qJD(1), 0, qJD(1) * t122, t124 * t138, t118 * r_i_i_C(1) - t119 * r_i_i_C(2); t122 * t137 + t119 * r_i_i_C(1) + t118 * r_i_i_C(2) - t124 * qJD(3) + (cos(qJ(1)) * pkin(1) + t135 * t122 + t142 * t124) * qJD(1), 0, -qJD(1) * t124, t122 * t138, -t116 * r_i_i_C(1) + t117 * r_i_i_C(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end