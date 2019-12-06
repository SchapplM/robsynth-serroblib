% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRRP2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:02:09
	% EndTime: 2019-12-05 18:02:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:02:09
	% EndTime: 2019-12-05 18:02:09
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
	% StartTime: 2019-12-05 18:02:09
	% EndTime: 2019-12-05 18:02:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->2), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(8);
	t13 = qJD(1) * cos(t12);
	t10 = qJD(1) * sin(t12);
	t1 = [0, 0, 0, 0, 0; t10, 0, 0, 0, 0; -t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; t13, 0, 0, 0, 0; t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:02:09
	% EndTime: 2019-12-05 18:02:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (26->4), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t27 = qJD(1) + qJD(3);
	t26 = qJ(1) + pkin(8) + qJ(3);
	t24 = t27 * cos(t26);
	t23 = t27 * sin(t26);
	t1 = [0, 0, 0, 0, 0; t23, 0, t23, 0, 0; -t24, 0, -t24, 0, 0; 0, 0, 0, 0, 0; t24, 0, t24, 0, 0; t23, 0, t23, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:02:09
	% EndTime: 2019-12-05 18:02:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (90->16), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t116 = qJ(1) + pkin(8) + qJ(3);
	t114 = sin(t116);
	t117 = qJD(1) + qJD(3);
	t125 = t117 * t114;
	t115 = cos(t116);
	t124 = t117 * t115;
	t118 = sin(qJ(4));
	t123 = t117 * t118;
	t119 = cos(qJ(4));
	t122 = t117 * t119;
	t121 = qJD(4) * t118;
	t120 = qJD(4) * t119;
	t111 = -t114 * t121 + t115 * t122;
	t110 = t114 * t120 + t115 * t123;
	t109 = t114 * t122 + t115 * t121;
	t108 = t114 * t123 - t115 * t120;
	t1 = [0, 0, 0, -t121, 0; t109, 0, t109, t110, 0; -t111, 0, -t111, t108, 0; 0, 0, 0, -t120, 0; -t108, 0, -t108, t111, 0; t110, 0, t110, t109, 0; 0, 0, 0, 0, 0; -t124, 0, -t124, 0, 0; -t125, 0, -t125, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:02:09
	% EndTime: 2019-12-05 18:02:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (90->16), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t139 = qJ(1) + pkin(8) + qJ(3);
	t137 = sin(t139);
	t140 = qJD(1) + qJD(3);
	t148 = t140 * t137;
	t138 = cos(t139);
	t147 = t140 * t138;
	t141 = sin(qJ(4));
	t146 = t140 * t141;
	t142 = cos(qJ(4));
	t145 = t140 * t142;
	t144 = qJD(4) * t141;
	t143 = qJD(4) * t142;
	t134 = -t137 * t144 + t138 * t145;
	t133 = t137 * t143 + t138 * t146;
	t132 = t137 * t145 + t138 * t144;
	t131 = t137 * t146 - t138 * t143;
	t1 = [0, 0, 0, -t144, 0; t132, 0, t132, t133, 0; -t134, 0, -t134, t131, 0; 0, 0, 0, -t143, 0; -t131, 0, -t131, t134, 0; t133, 0, t133, t132, 0; 0, 0, 0, 0, 0; -t147, 0, -t147, 0, 0; -t148, 0, -t148, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end