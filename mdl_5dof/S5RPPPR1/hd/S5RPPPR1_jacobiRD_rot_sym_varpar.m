% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPPR1
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
%   Siehe auch: S5RPPPR1_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:23
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:23
	% EndTime: 2022-01-20 09:13:23
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
	% StartTime: 2022-01-20 09:13:22
	% EndTime: 2022-01-20 09:13:22
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(7);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0; -t37, 0, 0, 0, 0; 0, 0, 0, 0, 0; t37, 0, 0, 0, 0; -t36, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:23
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (9->4), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t26 = qJD(1) * sin(pkin(8));
	t25 = qJD(1) * cos(pkin(8));
	t22 = qJ(1) + pkin(7);
	t21 = cos(t22);
	t20 = sin(t22);
	t1 = [-t21 * t25, 0, 0, 0, 0; -t20 * t25, 0, 0, 0, 0; 0, 0, 0, 0, 0; t21 * t26, 0, 0, 0, 0; t20 * t26, 0, 0, 0, 0; 0, 0, 0, 0, 0; -qJD(1) * t20, 0, 0, 0, 0; qJD(1) * t21, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:23
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (16->7), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->10)
	t125 = sin(pkin(9));
	t128 = cos(pkin(8));
	t131 = t125 * t128;
	t127 = cos(pkin(9));
	t130 = t127 * t128;
	t129 = qJD(1) * sin(pkin(8));
	t124 = qJ(1) + pkin(7);
	t123 = cos(t124);
	t122 = sin(t124);
	t1 = [(-t122 * t125 - t123 * t130) * qJD(1), 0, 0, 0, 0; (-t122 * t130 + t123 * t125) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; (-t122 * t127 + t123 * t131) * qJD(1), 0, 0, 0, 0; (t122 * t131 + t123 * t127) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; -t123 * t129, 0, 0, 0, 0; -t122 * t129, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:23
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (95->15), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->21)
	t178 = qJ(1) + pkin(7);
	t174 = sin(t178);
	t180 = cos(pkin(8));
	t188 = t174 * t180;
	t176 = cos(t178);
	t187 = t176 * t180;
	t179 = sin(pkin(8));
	t186 = qJD(1) * t179;
	t185 = qJD(5) * t179;
	t177 = pkin(9) + qJ(5);
	t173 = sin(t177);
	t175 = cos(t177);
	t184 = -t173 * t174 - t175 * t187;
	t183 = -t173 * t176 + t175 * t188;
	t182 = t173 * t187 - t174 * t175;
	t181 = t173 * t188 + t175 * t176;
	t172 = t184 * qJD(1) + t181 * qJD(5);
	t171 = t182 * qJD(1) + t183 * qJD(5);
	t170 = t183 * qJD(1) + t182 * qJD(5);
	t169 = t181 * qJD(1) + t184 * qJD(5);
	t1 = [t172, 0, 0, 0, t169; -t170, 0, 0, 0, -t171; 0, 0, 0, 0, -t175 * t185; t171, 0, 0, 0, t170; t169, 0, 0, 0, t172; 0, 0, 0, 0, t173 * t185; -t176 * t186, 0, 0, 0, 0; -t174 * t186, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end