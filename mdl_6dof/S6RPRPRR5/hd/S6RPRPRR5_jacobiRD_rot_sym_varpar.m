% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:53
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
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
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t22 = qJD(1) * sin(qJ(1));
	t21 = qJD(1) * cos(qJ(1));
	t18 = cos(pkin(10));
	t17 = sin(pkin(10));
	t1 = [-t18 * t21, 0, 0, 0, 0, 0; -t18 * t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t17 * t21, 0, 0, 0, 0, 0; t17 * t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t22, 0, 0, 0, 0, 0; t21, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t43 = sin(qJ(1));
	t48 = qJD(1) * t43;
	t44 = cos(qJ(1));
	t47 = qJD(1) * t44;
	t46 = qJD(3) * t43;
	t45 = qJD(3) * t44;
	t42 = pkin(10) + qJ(3);
	t41 = cos(t42);
	t40 = sin(t42);
	t39 = t40 * t46 - t41 * t47;
	t38 = t40 * t47 + t41 * t46;
	t37 = t40 * t45 + t41 * t48;
	t36 = t40 * t48 - t41 * t45;
	t1 = [t39, 0, t36, 0, 0, 0; -t37, 0, -t38, 0, 0, 0; 0, 0, -qJD(3) * t40, 0, 0, 0; t38, 0, t37, 0, 0, 0; t36, 0, t39, 0, 0, 0; 0, 0, -qJD(3) * t41, 0, 0, 0; -t48, 0, 0, 0, 0, 0; t47, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:11
	% EndTime: 2019-10-10 00:53:11
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (28->9), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t166 = sin(qJ(1));
	t171 = qJD(1) * t166;
	t167 = cos(qJ(1));
	t170 = qJD(1) * t167;
	t169 = qJD(3) * t166;
	t168 = qJD(3) * t167;
	t165 = pkin(10) + qJ(3);
	t164 = cos(t165);
	t163 = sin(t165);
	t162 = -t163 * t169 + t164 * t170;
	t161 = -t163 * t170 - t164 * t169;
	t160 = -t163 * t168 - t164 * t171;
	t159 = t163 * t171 - t164 * t168;
	t1 = [-t162, 0, t159, 0, 0, 0; t160, 0, t161, 0, 0, 0; 0, 0, -qJD(3) * t163, 0, 0, 0; -t171, 0, 0, 0, 0, 0; t170, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t161, 0, t160, 0, 0, 0; -t159, 0, t162, 0, 0, 0; 0, 0, qJD(3) * t164, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (161->17), mult. (250->18), div. (0->0), fcn. (250->6), ass. (0->19)
	t132 = qJD(3) - qJD(5);
	t116 = pkin(10) + qJ(3);
	t114 = sin(t116);
	t115 = cos(t116);
	t117 = sin(qJ(5));
	t119 = cos(qJ(5));
	t133 = -t114 * t119 + t115 * t117;
	t134 = t132 * t133;
	t123 = t114 * t117 + t115 * t119;
	t109 = t132 * t123;
	t122 = qJD(1) * t133;
	t121 = qJD(1) * t123;
	t120 = cos(qJ(1));
	t118 = sin(qJ(1));
	t107 = t118 * t134 + t120 * t121;
	t106 = -t118 * t109 + t120 * t122;
	t105 = t118 * t121 - t120 * t134;
	t104 = t109 * t120 + t118 * t122;
	t1 = [-t107, 0, -t104, 0, t104, 0; -t105, 0, t106, 0, -t106, 0; 0, 0, t134, 0, -t134, 0; t106, 0, -t105, 0, t105, 0; t104, 0, t107, 0, -t107, 0; 0, 0, t109, 0, -t109, 0; qJD(1) * t118, 0, 0, 0, 0, 0; -qJD(1) * t120, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:13
	% EndTime: 2019-10-10 00:53:13
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (402->47), mult. (634->69), div. (0->0), fcn. (656->8), ass. (0->43)
	t445 = pkin(10) + qJ(3);
	t443 = sin(t445);
	t444 = cos(t445);
	t447 = sin(qJ(5));
	t450 = cos(qJ(5));
	t484 = -t443 * t450 + t444 * t447;
	t476 = qJD(3) * t447;
	t487 = -t484 * qJD(5) + t444 * t476;
	t424 = t443 * t447 + t444 * t450;
	t451 = cos(qJ(1));
	t423 = t424 * t451;
	t475 = qJD(3) * t450;
	t419 = -t443 * t475 + t487;
	t486 = t424 * qJD(5) - t443 * t476;
	t485 = t484 * qJD(1);
	t448 = sin(qJ(1));
	t456 = qJD(1) * t424;
	t467 = t448 * t475;
	t417 = -t443 * t467 + t487 * t448 + t451 * t456;
	t421 = t424 * t448;
	t446 = sin(qJ(6));
	t449 = cos(qJ(6));
	t478 = qJD(1) * t448;
	t483 = (t421 * t449 + t446 * t451) * qJD(6) + t417 * t446 + t449 * t478;
	t477 = qJD(1) * t451;
	t473 = qJD(6) * t446;
	t472 = qJD(6) * t449;
	t415 = t485 * t448 + (qJD(3) - qJD(5)) * t423;
	t422 = t484 * t451;
	t462 = -t415 * t446 + t422 * t472;
	t461 = t415 * t449 + t422 * t473;
	t416 = -t444 * t467 + t486 * t448 + t485 * t451;
	t420 = t484 * t448;
	t460 = t416 * t446 + t420 * t472;
	t459 = -t416 * t449 + t420 * t473;
	t458 = t419 * t446 + t424 * t472;
	t457 = -t419 * t449 + t424 * t473;
	t452 = -t417 * t449 + t446 * t478 + (t421 * t446 - t449 * t451) * qJD(6);
	t418 = -t444 * t475 + t486;
	t414 = -t419 * t451 + t448 * t456;
	t413 = -t446 * t477 - t414 * t449 + (-t423 * t446 - t448 * t449) * qJD(6);
	t412 = -t449 * t477 + t414 * t446 + (-t423 * t449 + t446 * t448) * qJD(6);
	t1 = [t452, 0, -t461, 0, t461, t412; t413, 0, -t459, 0, t459, -t483; 0, 0, -t457, 0, t457, t418 * t446 + t472 * t484; t483, 0, -t462, 0, t462, -t413; t412, 0, -t460, 0, t460, t452; 0, 0, -t458, 0, t458, t418 * t449 - t473 * t484; -t416, 0, t414, 0, -t414, 0; -t415, 0, -t417, 0, t417, 0; 0, 0, t418, 0, -t418, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end