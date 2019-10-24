% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:47
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:47:13
	% EndTime: 2019-10-24 10:47:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:47:13
	% EndTime: 2019-10-24 10:47:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = qJD(1) * cos(qJ(1));
	t7 = qJD(1) * sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t7, 0, 0, 0, 0; -t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9, 0, 0, 0, 0; t7, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:47:13
	% EndTime: 2019-10-24 10:47:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (18->4), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t26 = qJ(1) + qJ(2);
	t25 = qJD(1) + qJD(2);
	t23 = t25 * cos(t26);
	t22 = t25 * sin(t26);
	t1 = [0, 0, 0, 0, 0; t22, t22, 0, 0, 0; -t23, -t23, 0, 0, 0; 0, 0, 0, 0, 0; t23, t23, 0, 0, 0; t22, t22, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:47:13
	% EndTime: 2019-10-24 10:47:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (26->4), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t31 = qJD(1) + qJD(2);
	t30 = qJ(1) + qJ(2) + pkin(8);
	t28 = t31 * cos(t30);
	t27 = t31 * sin(t30);
	t1 = [0, 0, 0, 0, 0; t27, t27, 0, 0, 0; -t28, -t28, 0, 0, 0; 0, 0, 0, 0, 0; t28, t28, 0, 0, 0; t27, t27, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:47:13
	% EndTime: 2019-10-24 10:47:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (44->10), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t80 = qJ(1) + qJ(2) + pkin(8);
	t78 = sin(t80);
	t81 = qJD(1) + qJD(2);
	t89 = t81 * t78;
	t79 = cos(t80);
	t88 = t81 * t79;
	t87 = t81 * sin(pkin(9));
	t86 = t81 * cos(pkin(9));
	t85 = t78 * t87;
	t84 = t79 * t86;
	t77 = t79 * t87;
	t76 = t78 * t86;
	t1 = [0, 0, 0, 0, 0; t76, t76, 0, 0, 0; -t84, -t84, 0, 0, 0; 0, 0, 0, 0, 0; -t85, -t85, 0, 0, 0; t77, t77, 0, 0, 0; 0, 0, 0, 0, 0; -t88, -t88, 0, 0, 0; -t89, -t89, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:47:13
	% EndTime: 2019-10-24 10:47:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (116->17), mult. (54->12), div. (0->0), fcn. (54->4), ass. (0->16)
	t136 = qJ(1) + qJ(2) + pkin(8);
	t132 = sin(t136);
	t138 = qJD(1) + qJD(2);
	t142 = t138 * t132;
	t133 = cos(t136);
	t141 = t138 * t133;
	t137 = pkin(9) + qJ(5);
	t134 = sin(t137);
	t140 = qJD(5) * t134;
	t135 = cos(t137);
	t139 = qJD(5) * t135;
	t129 = -t132 * t140 + t135 * t141;
	t128 = t132 * t139 + t134 * t141;
	t127 = t133 * t140 + t135 * t142;
	t126 = -t133 * t139 + t134 * t142;
	t1 = [0, 0, 0, 0, -t140; t127, t127, 0, 0, t128; -t129, -t129, 0, 0, t126; 0, 0, 0, 0, -t139; -t126, -t126, 0, 0, t129; t128, t128, 0, 0, t127; 0, 0, 0, 0, 0; -t141, -t141, 0, 0, 0; -t142, -t142, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end