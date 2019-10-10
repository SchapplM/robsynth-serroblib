% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:15
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRPR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.07s
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
	% StartTime: 2019-10-10 10:15:25
	% EndTime: 2019-10-10 10:15:25
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
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.10s
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
	t1 = [-t98, -t95, 0, t95, 0, 0; -t96, t97, 0, -t97, 0, 0; 0, -t124, 0, t124, 0, 0; t97, -t96, 0, t96, 0, 0; t95, t98, 0, -t98, 0, 0; 0, t100, 0, -t100, 0, 0; qJD(1) * t107, 0, 0, 0, 0, 0; -qJD(1) * t110, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (161->17), mult. (250->18), div. (0->0), fcn. (250->6), ass. (0->19)
	t141 = qJD(2) - qJD(4);
	t125 = qJ(4) + pkin(10);
	t123 = sin(t125);
	t124 = cos(t125);
	t126 = sin(qJ(2));
	t128 = cos(qJ(2));
	t142 = -t123 * t128 + t124 * t126;
	t143 = t141 * t142;
	t132 = t123 * t126 + t124 * t128;
	t118 = t141 * t132;
	t131 = qJD(1) * t142;
	t130 = qJD(1) * t132;
	t129 = cos(qJ(1));
	t127 = sin(qJ(1));
	t116 = -t127 * t143 + t129 * t130;
	t115 = -t127 * t118 - t129 * t131;
	t114 = t127 * t130 + t129 * t143;
	t113 = t118 * t129 - t127 * t131;
	t1 = [-t116, -t113, 0, t113, 0, 0; -t114, t115, 0, -t115, 0, 0; 0, -t143, 0, t143, 0, 0; t115, -t114, 0, t114, 0, 0; t113, t116, 0, -t116, 0, 0; 0, t118, 0, -t118, 0, 0; qJD(1) * t127, 0, 0, 0, 0, 0; -qJD(1) * t129, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:27
	% EndTime: 2019-10-10 10:15:27
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (402->47), mult. (634->69), div. (0->0), fcn. (656->8), ass. (0->43)
	t458 = qJ(4) + pkin(10);
	t456 = sin(t458);
	t457 = cos(t458);
	t460 = sin(qJ(2));
	t463 = cos(qJ(2));
	t497 = t463 * t456 - t460 * t457;
	t488 = qJD(2) * t463;
	t500 = -t497 * qJD(4) + t456 * t488;
	t437 = t460 * t456 + t463 * t457;
	t464 = cos(qJ(1));
	t436 = t437 * t464;
	t489 = qJD(2) * t460;
	t432 = -t457 * t489 + t500;
	t499 = t437 * qJD(4) - t457 * t488;
	t498 = t497 * qJD(1);
	t461 = sin(qJ(1));
	t469 = qJD(1) * t437;
	t482 = t461 * t489;
	t430 = -t457 * t482 + t500 * t461 + t464 * t469;
	t434 = t437 * t461;
	t459 = sin(qJ(6));
	t462 = cos(qJ(6));
	t491 = qJD(1) * t461;
	t496 = (t434 * t462 + t459 * t464) * qJD(6) + t430 * t459 + t462 * t491;
	t490 = qJD(1) * t464;
	t486 = qJD(6) * t459;
	t485 = qJD(6) * t462;
	t428 = t498 * t461 + (qJD(2) - qJD(4)) * t436;
	t435 = t497 * t464;
	t475 = -t428 * t459 + t435 * t485;
	t474 = t428 * t462 + t435 * t486;
	t429 = -t456 * t482 + t499 * t461 + t498 * t464;
	t433 = t497 * t461;
	t473 = t429 * t459 + t433 * t485;
	t472 = -t429 * t462 + t433 * t486;
	t471 = t432 * t459 + t437 * t485;
	t470 = -t432 * t462 + t437 * t486;
	t465 = -t430 * t462 + t459 * t491 + (t434 * t459 - t462 * t464) * qJD(6);
	t431 = -t456 * t489 + t499;
	t427 = -t432 * t464 + t461 * t469;
	t426 = -t459 * t490 - t427 * t462 + (-t436 * t459 - t461 * t462) * qJD(6);
	t425 = -t462 * t490 + t427 * t459 + (-t436 * t462 + t459 * t461) * qJD(6);
	t1 = [t465, -t474, 0, t474, 0, t425; t426, -t472, 0, t472, 0, -t496; 0, -t470, 0, t470, 0, t431 * t459 + t485 * t497; t496, -t475, 0, t475, 0, -t426; t425, -t473, 0, t473, 0, t465; 0, -t471, 0, t471, 0, t431 * t462 - t486 * t497; -t429, t427, 0, -t427, 0, 0; -t428, -t430, 0, t430, 0, 0; 0, t431, 0, -t431, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end