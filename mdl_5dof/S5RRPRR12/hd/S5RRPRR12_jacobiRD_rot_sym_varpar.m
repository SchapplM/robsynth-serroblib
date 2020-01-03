% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPRR12_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR12_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:31:15
	% EndTime: 2019-12-31 20:31:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:31:15
	% EndTime: 2019-12-31 20:31:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:31:15
	% EndTime: 2019-12-31 20:31:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t34 = sin(qJ(1));
	t41 = qJD(1) * t34;
	t36 = cos(qJ(1));
	t40 = qJD(1) * t36;
	t33 = sin(qJ(2));
	t39 = qJD(2) * t33;
	t35 = cos(qJ(2));
	t38 = qJD(2) * t35;
	t37 = qJD(2) * t36;
	t32 = t34 * t39 - t35 * t40;
	t31 = t33 * t40 + t34 * t38;
	t30 = t33 * t37 + t35 * t41;
	t29 = t33 * t41 - t35 * t37;
	t1 = [t32, t29, 0, 0, 0; -t30, -t31, 0, 0, 0; 0, -t39, 0, 0, 0; t31, t30, 0, 0, 0; t29, t32, 0, 0, 0; 0, -t38, 0, 0, 0; -t41, 0, 0, 0, 0; t40, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:31:15
	% EndTime: 2019-12-31 20:31:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->8), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t156 = sin(qJ(1));
	t163 = qJD(1) * t156;
	t158 = cos(qJ(1));
	t162 = qJD(1) * t158;
	t155 = sin(qJ(2));
	t161 = qJD(2) * t155;
	t157 = cos(qJ(2));
	t160 = qJD(2) * t157;
	t159 = qJD(2) * t158;
	t154 = -t156 * t161 + t157 * t162;
	t153 = -t155 * t162 - t156 * t160;
	t152 = -t155 * t159 - t157 * t163;
	t151 = t155 * t163 - t157 * t159;
	t1 = [-t154, t151, 0, 0, 0; t152, t153, 0, 0, 0; 0, -t161, 0, 0, 0; -t163, 0, 0, 0, 0; t162, 0, 0, 0, 0; 0, 0, 0, 0, 0; t153, t152, 0, 0, 0; -t151, t154, 0, 0, 0; 0, t160, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:31:15
	% EndTime: 2019-12-31 20:31:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (73->16), mult. (250->18), div. (0->0), fcn. (250->6), ass. (0->18)
	t122 = qJD(2) - qJD(4);
	t105 = sin(qJ(4));
	t106 = sin(qJ(2));
	t108 = cos(qJ(4));
	t109 = cos(qJ(2));
	t123 = -t105 * t109 + t106 * t108;
	t124 = t122 * t123;
	t113 = t105 * t106 + t108 * t109;
	t100 = t122 * t113;
	t112 = qJD(1) * t123;
	t111 = qJD(1) * t113;
	t110 = cos(qJ(1));
	t107 = sin(qJ(1));
	t98 = -t107 * t124 + t110 * t111;
	t97 = -t107 * t100 - t110 * t112;
	t96 = t107 * t111 + t110 * t124;
	t95 = t100 * t110 - t107 * t112;
	t1 = [-t98, -t95, 0, t95, 0; -t96, t97, 0, -t97, 0; 0, -t124, 0, t124, 0; t97, -t96, 0, t96, 0; t95, t98, 0, -t98, 0; 0, t100, 0, -t100, 0; qJD(1) * t107, 0, 0, 0, 0; -qJD(1) * t110, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:31:17
	% EndTime: 2019-12-31 20:31:17
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (194->44), mult. (634->66), div. (0->0), fcn. (656->8), ass. (0->39)
	t435 = sin(qJ(4));
	t436 = sin(qJ(2));
	t439 = cos(qJ(4));
	t440 = cos(qJ(2));
	t475 = t440 * t435 - t436 * t439;
	t415 = t436 * t435 + t440 * t439;
	t441 = cos(qJ(1));
	t414 = t415 * t441;
	t466 = qJD(2) * t440;
	t467 = qJD(2) * t436;
	t410 = -t475 * qJD(4) + t435 * t466 - t439 * t467;
	t476 = t475 * qJD(1);
	t409 = t415 * qJD(4) - t435 * t467 - t439 * t466;
	t437 = sin(qJ(1));
	t446 = qJD(1) * t415;
	t408 = t410 * t437 + t441 * t446;
	t412 = t415 * t437;
	t434 = sin(qJ(5));
	t438 = cos(qJ(5));
	t469 = qJD(1) * t437;
	t474 = (t412 * t438 + t434 * t441) * qJD(5) + t408 * t434 + t438 * t469;
	t468 = qJD(1) * t441;
	t464 = qJD(5) * t434;
	t463 = qJD(5) * t438;
	t406 = t476 * t437 + (qJD(2) - qJD(4)) * t414;
	t413 = t475 * t441;
	t452 = -t406 * t434 + t413 * t463;
	t451 = t406 * t438 + t413 * t464;
	t407 = t409 * t437 + t476 * t441;
	t411 = t475 * t437;
	t450 = t407 * t434 + t411 * t463;
	t449 = -t407 * t438 + t411 * t464;
	t448 = t410 * t434 + t415 * t463;
	t447 = -t410 * t438 + t415 * t464;
	t442 = -t408 * t438 + t434 * t469 + (t412 * t434 - t438 * t441) * qJD(5);
	t405 = -t410 * t441 + t437 * t446;
	t404 = -t434 * t468 - t405 * t438 + (-t414 * t434 - t437 * t438) * qJD(5);
	t403 = -t438 * t468 + t405 * t434 + (-t414 * t438 + t434 * t437) * qJD(5);
	t1 = [t442, -t451, 0, t451, t403; t404, -t449, 0, t449, -t474; 0, -t447, 0, t447, t409 * t434 + t463 * t475; t474, -t452, 0, t452, -t404; t403, -t450, 0, t450, t442; 0, -t448, 0, t448, t409 * t438 - t464 * t475; -t407, t405, 0, -t405, 0; -t406, -t408, 0, t408, 0; 0, t409, 0, -t409, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end