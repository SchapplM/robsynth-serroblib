% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPPRR6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.06s
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
	t1 = [t32, t29, 0, 0, 0, 0; -t30, -t31, 0, 0, 0, 0; 0, -t39, 0, 0, 0, 0; t31, t30, 0, 0, 0, 0; t29, t32, 0, 0, 0, 0; 0, -t38, 0, 0, 0, 0; -t41, 0, 0, 0, 0, 0; t40, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:55
	% EndTime: 2019-10-10 09:44:55
	% DurationCPUTime: 0.04s
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
	t1 = [-t154, t151, 0, 0, 0, 0; t152, t153, 0, 0, 0, 0; 0, -t161, 0, 0, 0, 0; -t163, 0, 0, 0, 0, 0; t162, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t153, t152, 0, 0, 0, 0; -t151, t154, 0, 0, 0, 0; 0, t160, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (27->9), mult. (106->18), div. (0->0), fcn. (106->6), ass. (0->18)
	t67 = cos(qJ(1));
	t74 = qJD(1) * t67;
	t62 = sin(pkin(10));
	t63 = cos(pkin(10));
	t64 = sin(qJ(2));
	t66 = cos(qJ(2));
	t73 = t62 * t66 - t63 * t64;
	t72 = t62 * t64 + t63 * t66;
	t65 = sin(qJ(1));
	t71 = t72 * t65;
	t70 = qJD(1) * t73;
	t69 = t73 * qJD(2);
	t68 = t72 * qJD(2);
	t61 = t65 * t69 + t72 * t74;
	t60 = -qJD(2) * t71 + t67 * t70;
	t59 = -qJD(1) * t71 + t67 * t69;
	t58 = t65 * t70 + t67 * t68;
	t1 = [-t61, -t58, 0, 0, 0, 0; t59, t60, 0, 0, 0, 0; 0, t69, 0, 0, 0, 0; t60, t59, 0, 0, 0, 0; t58, t61, 0, 0, 0, 0; 0, t68, 0, 0, 0, 0; qJD(1) * t65, 0, 0, 0, 0, 0; -t74, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (161->17), mult. (250->18), div. (0->0), fcn. (250->6), ass. (0->19)
	t132 = qJD(2) - qJD(5);
	t116 = pkin(10) + qJ(5);
	t114 = sin(t116);
	t115 = cos(t116);
	t117 = sin(qJ(2));
	t119 = cos(qJ(2));
	t133 = -t114 * t119 + t115 * t117;
	t134 = t132 * t133;
	t123 = t114 * t117 + t115 * t119;
	t109 = t132 * t123;
	t122 = qJD(1) * t133;
	t121 = qJD(1) * t123;
	t120 = cos(qJ(1));
	t118 = sin(qJ(1));
	t107 = -t118 * t134 + t120 * t121;
	t106 = -t118 * t109 - t120 * t122;
	t105 = t118 * t121 + t120 * t134;
	t104 = t109 * t120 - t118 * t122;
	t1 = [-t107, -t104, 0, 0, t104, 0; -t105, t106, 0, 0, -t106, 0; 0, -t134, 0, 0, t134, 0; t106, -t105, 0, 0, t105, 0; t104, t107, 0, 0, -t107, 0; 0, t109, 0, 0, -t109, 0; qJD(1) * t118, 0, 0, 0, 0, 0; -qJD(1) * t120, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:57
	% EndTime: 2019-10-10 09:44:57
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (402->47), mult. (634->69), div. (0->0), fcn. (656->8), ass. (0->43)
	t449 = pkin(10) + qJ(5);
	t447 = sin(t449);
	t448 = cos(t449);
	t451 = sin(qJ(2));
	t454 = cos(qJ(2));
	t488 = t454 * t447 - t451 * t448;
	t479 = qJD(2) * t454;
	t491 = -t488 * qJD(5) + t447 * t479;
	t428 = t451 * t447 + t454 * t448;
	t455 = cos(qJ(1));
	t427 = t428 * t455;
	t480 = qJD(2) * t451;
	t423 = -t448 * t480 + t491;
	t490 = t428 * qJD(5) - t448 * t479;
	t489 = t488 * qJD(1);
	t452 = sin(qJ(1));
	t460 = qJD(1) * t428;
	t473 = t452 * t480;
	t421 = -t448 * t473 + t491 * t452 + t455 * t460;
	t425 = t428 * t452;
	t450 = sin(qJ(6));
	t453 = cos(qJ(6));
	t482 = qJD(1) * t452;
	t487 = (t425 * t453 + t450 * t455) * qJD(6) + t421 * t450 + t453 * t482;
	t481 = qJD(1) * t455;
	t477 = qJD(6) * t450;
	t476 = qJD(6) * t453;
	t419 = t489 * t452 + (qJD(2) - qJD(5)) * t427;
	t426 = t488 * t455;
	t466 = -t419 * t450 + t426 * t476;
	t465 = t419 * t453 + t426 * t477;
	t420 = -t447 * t473 + t490 * t452 + t489 * t455;
	t424 = t488 * t452;
	t464 = t420 * t450 + t424 * t476;
	t463 = -t420 * t453 + t424 * t477;
	t462 = t423 * t450 + t428 * t476;
	t461 = -t423 * t453 + t428 * t477;
	t456 = -t421 * t453 + t450 * t482 + (t425 * t450 - t453 * t455) * qJD(6);
	t422 = -t447 * t480 + t490;
	t418 = -t423 * t455 + t452 * t460;
	t417 = -t450 * t481 - t418 * t453 + (-t427 * t450 - t452 * t453) * qJD(6);
	t416 = -t453 * t481 + t418 * t450 + (-t427 * t453 + t450 * t452) * qJD(6);
	t1 = [t456, -t465, 0, 0, t465, t416; t417, -t463, 0, 0, t463, -t487; 0, -t461, 0, 0, t461, t422 * t450 + t476 * t488; t487, -t466, 0, 0, t466, -t417; t416, -t464, 0, 0, t464, t456; 0, -t462, 0, 0, t462, t422 * t453 - t477 * t488; -t420, t418, 0, 0, -t418, 0; -t419, -t421, 0, 0, t421, 0; 0, t422, 0, 0, -t422, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end