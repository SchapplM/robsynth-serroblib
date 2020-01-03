% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPRP8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:04:51
	% EndTime: 2019-12-31 20:04:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:04:51
	% EndTime: 2019-12-31 20:04:51
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
	% StartTime: 2019-12-31 20:04:51
	% EndTime: 2019-12-31 20:04:51
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
	% StartTime: 2019-12-31 20:04:51
	% EndTime: 2019-12-31 20:04:51
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
	% StartTime: 2019-12-31 20:04:51
	% EndTime: 2019-12-31 20:04:51
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
	% StartTime: 2019-12-31 20:04:51
	% EndTime: 2019-12-31 20:04:51
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (73->16), mult. (250->18), div. (0->0), fcn. (250->6), ass. (0->18)
	t135 = qJD(2) - qJD(4);
	t118 = sin(qJ(4));
	t119 = sin(qJ(2));
	t121 = cos(qJ(4));
	t122 = cos(qJ(2));
	t136 = -t118 * t122 + t119 * t121;
	t137 = t135 * t136;
	t126 = t118 * t119 + t121 * t122;
	t113 = t135 * t126;
	t125 = qJD(1) * t136;
	t124 = qJD(1) * t126;
	t123 = cos(qJ(1));
	t120 = sin(qJ(1));
	t111 = -t120 * t137 + t123 * t124;
	t110 = -t113 * t120 - t123 * t125;
	t109 = t120 * t124 + t123 * t137;
	t108 = t113 * t123 - t120 * t125;
	t1 = [-t111, -t108, 0, t108, 0; -t109, t110, 0, -t110, 0; 0, -t137, 0, t137, 0; t110, -t109, 0, t109, 0; t108, t111, 0, -t111, 0; 0, t113, 0, -t113, 0; qJD(1) * t120, 0, 0, 0, 0; -qJD(1) * t123, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end