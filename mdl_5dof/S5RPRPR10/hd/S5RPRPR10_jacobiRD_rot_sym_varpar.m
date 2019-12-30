% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRPR10_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR10_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:57:18
	% EndTime: 2019-12-29 16:57:18
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:57:23
	% EndTime: 2019-12-29 16:57:23
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:57:18
	% EndTime: 2019-12-29 16:57:18
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t14 = qJD(1) * sin(qJ(1));
	t13 = qJD(1) * cos(qJ(1));
	t1 = [-t13, 0, 0, 0, 0; -t14, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t14, 0, 0, 0, 0; t13, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:57:23
	% EndTime: 2019-12-29 16:57:23
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (24->10), mult. (64->12), div. (0->0), fcn. (64->4), ass. (0->9)
	t95 = sin(qJ(1));
	t99 = qJD(3) * t95;
	t97 = cos(qJ(1));
	t98 = qJD(3) * t97;
	t96 = cos(qJ(3));
	t94 = sin(qJ(3));
	t89 = -t94 * t98 + t96 * t99 + (t94 * t97 - t95 * t96) * qJD(1);
	t88 = -t94 * t99 - t96 * t98 + (t94 * t95 + t96 * t97) * qJD(1);
	t1 = [-t88, 0, t88, 0, 0; t89, 0, -t89, 0, 0; 0, 0, 0, 0, 0; t89, 0, -t89, 0, 0; t88, 0, -t88, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:57:18
	% EndTime: 2019-12-29 16:57:18
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (56->11), mult. (64->12), div. (0->0), fcn. (64->4), ass. (0->10)
	t115 = sin(qJ(1));
	t118 = qJD(3) * t115;
	t116 = cos(qJ(1));
	t117 = qJD(3) * t116;
	t114 = qJ(3) + pkin(8);
	t113 = cos(t114);
	t112 = sin(t114);
	t107 = -t112 * t117 + t113 * t118 + (t112 * t116 - t113 * t115) * qJD(1);
	t106 = -t112 * t118 - t113 * t117 + (t112 * t115 + t113 * t116) * qJD(1);
	t1 = [-t106, 0, t106, 0, 0; t107, 0, -t107, 0, 0; 0, 0, 0, 0, 0; t107, 0, -t107, 0, 0; t106, 0, -t106, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:57:13
	% EndTime: 2019-12-29 16:57:13
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (160->16), mult. (190->16), div. (0->0), fcn. (202->6), ass. (0->19)
	t126 = qJD(1) - qJD(3);
	t121 = qJ(3) + pkin(8);
	t119 = sin(t121);
	t120 = cos(t121);
	t124 = sin(qJ(1));
	t125 = cos(qJ(1));
	t104 = -t124 * t119 - t125 * t120;
	t105 = t125 * t119 - t124 * t120;
	t113 = sin(qJ(5));
	t123 = qJD(5) * t113;
	t114 = cos(qJ(5));
	t122 = qJD(5) * t114;
	t102 = t126 * t104;
	t98 = t102 * t113 + t105 * t122;
	t99 = -t102 * t114 + t105 * t123;
	t103 = t126 * t105;
	t100 = -t103 * t113 + t104 * t122;
	t101 = t103 * t114 + t104 * t123;
	t1 = [-t99, 0, t99, 0, t100; t101, 0, -t101, 0, t98; 0, 0, 0, 0, t123; -t98, 0, t98, 0, -t101; t100, 0, -t100, 0, -t99; 0, 0, 0, 0, t122; -t103, 0, t103, 0, 0; t102, 0, -t102, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end