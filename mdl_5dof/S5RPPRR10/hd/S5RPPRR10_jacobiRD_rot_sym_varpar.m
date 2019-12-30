% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRR10
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPRR10_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR10_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:24:03
	% EndTime: 2019-12-29 16:24:03
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:24:03
	% EndTime: 2019-12-29 16:24:03
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
	% StartTime: 2019-12-29 16:24:09
	% EndTime: 2019-12-29 16:24:09
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (3->3), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t22 = qJD(1) * sin(qJ(1));
	t21 = qJD(1) * cos(qJ(1));
	t18 = cos(pkin(8));
	t17 = sin(pkin(8));
	t1 = [-t18 * t21, 0, 0, 0, 0; -t18 * t22, 0, 0, 0, 0; 0, 0, 0, 0, 0; t17 * t21, 0, 0, 0, 0; t17 * t22, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t22, 0, 0, 0, 0; t21, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:24:04
	% EndTime: 2019-12-29 16:24:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t104 = qJD(1) * sin(qJ(1));
	t103 = qJD(1) * cos(qJ(1));
	t100 = cos(pkin(8));
	t99 = sin(pkin(8));
	t1 = [-t100 * t103, 0, 0, 0, 0; -t100 * t104, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t104, 0, 0, 0, 0; t103, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t99 * t103, 0, 0, 0, 0; -t99 * t104, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:24:03
	% EndTime: 2019-12-29 16:24:03
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (27->9), mult. (106->18), div. (0->0), fcn. (106->6), ass. (0->18)
	t57 = sin(qJ(1));
	t66 = qJD(1) * t57;
	t54 = sin(pkin(8));
	t55 = cos(pkin(8));
	t56 = sin(qJ(4));
	t58 = cos(qJ(4));
	t65 = -t54 * t58 + t55 * t56;
	t64 = t54 * t56 + t55 * t58;
	t59 = cos(qJ(1));
	t63 = t64 * t59;
	t62 = qJD(1) * t65;
	t61 = t65 * qJD(4);
	t60 = t64 * qJD(4);
	t53 = -qJD(1) * t63 + t57 * t61;
	t52 = t57 * t60 + t59 * t62;
	t51 = t59 * t61 + t64 * t66;
	t50 = -qJD(4) * t63 + t57 * t62;
	t1 = [t53, 0, 0, t50, 0; -t51, 0, 0, -t52, 0; 0, 0, 0, t61, 0; t52, 0, 0, t51, 0; t50, 0, 0, t53, 0; 0, 0, 0, t60, 0; t66, 0, 0, 0, 0; -qJD(1) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:24:04
	% EndTime: 2019-12-29 16:24:04
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (129->14), mult. (162->22), div. (0->0), fcn. (162->6), ass. (0->24)
	t106 = qJ(4) + qJ(5);
	t104 = cos(t106);
	t107 = sin(pkin(8));
	t119 = t104 * t107;
	t108 = cos(pkin(8));
	t118 = t104 * t108;
	t105 = qJD(4) + qJD(5);
	t109 = sin(qJ(1));
	t117 = t105 * t109;
	t116 = qJD(1) * t109;
	t110 = cos(qJ(1));
	t115 = qJD(1) * t110;
	t103 = sin(t106);
	t114 = t103 * t108 - t119;
	t113 = t103 * t107 + t118;
	t101 = t114 * t105;
	t112 = t114 * t109;
	t111 = t113 * t110;
	t100 = t113 * t105;
	t99 = -qJD(1) * t111 + t105 * t112;
	t98 = t117 * t118 - t115 * t119 + (t107 * t117 + t108 * t115) * t103;
	t97 = t110 * t101 + t113 * t116;
	t96 = qJD(1) * t112 - t105 * t111;
	t1 = [t99, 0, 0, t96, t96; -t97, 0, 0, -t98, -t98; 0, 0, 0, t101, t101; t98, 0, 0, t97, t97; t96, 0, 0, t99, t99; 0, 0, 0, t100, t100; t116, 0, 0, 0, 0; -t115, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end