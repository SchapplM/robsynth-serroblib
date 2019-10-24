% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:48
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPRP1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:48:01
	% EndTime: 2019-10-24 10:48:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:48:01
	% EndTime: 2019-10-24 10:48:01
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
	% StartTime: 2019-10-24 10:48:01
	% EndTime: 2019-10-24 10:48:01
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
	% StartTime: 2019-10-24 10:48:01
	% EndTime: 2019-10-24 10:48:01
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
	% StartTime: 2019-10-24 10:48:02
	% EndTime: 2019-10-24 10:48:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (90->16), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t121 = qJ(1) + qJ(2) + pkin(8);
	t119 = sin(t121);
	t122 = qJD(1) + qJD(2);
	t130 = t122 * t119;
	t120 = cos(t121);
	t129 = t122 * t120;
	t123 = sin(qJ(4));
	t128 = t122 * t123;
	t124 = cos(qJ(4));
	t127 = t122 * t124;
	t126 = qJD(4) * t123;
	t125 = qJD(4) * t124;
	t116 = -t119 * t126 + t120 * t127;
	t115 = t119 * t125 + t120 * t128;
	t114 = t119 * t127 + t120 * t126;
	t113 = t119 * t128 - t120 * t125;
	t1 = [0, 0, 0, -t126, 0; t114, t114, 0, t115, 0; -t116, -t116, 0, t113, 0; 0, 0, 0, -t125, 0; -t113, -t113, 0, t116, 0; t115, t115, 0, t114, 0; 0, 0, 0, 0, 0; -t129, -t129, 0, 0, 0; -t130, -t130, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:48:01
	% EndTime: 2019-10-24 10:48:02
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (90->16), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t147 = qJ(1) + qJ(2) + pkin(8);
	t145 = sin(t147);
	t148 = qJD(1) + qJD(2);
	t156 = t148 * t145;
	t146 = cos(t147);
	t155 = t148 * t146;
	t149 = sin(qJ(4));
	t154 = t148 * t149;
	t150 = cos(qJ(4));
	t153 = t148 * t150;
	t152 = qJD(4) * t149;
	t151 = qJD(4) * t150;
	t142 = -t145 * t152 + t146 * t153;
	t141 = t145 * t151 + t146 * t154;
	t140 = t145 * t153 + t146 * t152;
	t139 = t145 * t154 - t146 * t151;
	t1 = [0, 0, 0, -t152, 0; t140, t140, 0, t141, 0; -t142, -t142, 0, t139, 0; 0, 0, 0, -t151, 0; -t139, -t139, 0, t142, 0; t141, t141, 0, t140, 0; 0, 0, 0, 0, 0; -t155, -t155, 0, 0, 0; -t156, -t156, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end