% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPRP1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t10 = qJD(1) * sin(qJ(1));
	t9 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t10, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(7);
	t14 = qJD(1) * sin(t12);
	t13 = qJD(1) * cos(t12);
	t1 = [0, 0, 0, 0, 0; -t14, 0, 0, 0, 0; t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0; -t14, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->3), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t49 = qJD(1) * sin(pkin(8));
	t48 = qJD(1) * cos(pkin(8));
	t45 = qJ(1) + pkin(7);
	t44 = cos(t45);
	t43 = sin(t45);
	t1 = [0, 0, 0, 0, 0; -t43 * t48, 0, 0, 0, 0; t44 * t48, 0, 0, 0, 0; 0, 0, 0, 0, 0; t43 * t49, 0, 0, 0, 0; -t44 * t49, 0, 0, 0, 0; 0, 0, 0, 0, 0; qJD(1) * t44, 0, 0, 0, 0; qJD(1) * t43, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (60->13), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->20)
	t138 = cos(pkin(8));
	t139 = sin(qJ(4));
	t148 = t138 * t139;
	t140 = cos(qJ(4));
	t147 = t138 * t140;
	t137 = sin(pkin(8));
	t146 = qJD(1) * t137;
	t145 = qJD(4) * t137;
	t136 = qJ(1) + pkin(7);
	t134 = sin(t136);
	t135 = cos(t136);
	t144 = t134 * t139 + t135 * t147;
	t143 = t134 * t140 - t135 * t148;
	t142 = -t134 * t147 + t135 * t139;
	t141 = t134 * t148 + t135 * t140;
	t133 = t144 * qJD(1) - t141 * qJD(4);
	t132 = t143 * qJD(1) + t142 * qJD(4);
	t131 = t142 * qJD(1) + t143 * qJD(4);
	t130 = t141 * qJD(1) - t144 * qJD(4);
	t1 = [0, 0, 0, -t140 * t145, 0; t131, 0, 0, t132, 0; t133, 0, 0, -t130, 0; 0, 0, 0, t139 * t145, 0; t130, 0, 0, -t133, 0; t132, 0, 0, t131, 0; 0, 0, 0, 0, 0; -t134 * t146, 0, 0, 0, 0; t135 * t146, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (60->13), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->20)
	t148 = cos(pkin(8));
	t149 = sin(qJ(4));
	t158 = t148 * t149;
	t150 = cos(qJ(4));
	t157 = t148 * t150;
	t147 = sin(pkin(8));
	t156 = qJD(1) * t147;
	t155 = qJD(4) * t147;
	t146 = qJ(1) + pkin(7);
	t144 = sin(t146);
	t145 = cos(t146);
	t154 = t144 * t149 + t145 * t157;
	t153 = t144 * t150 - t145 * t158;
	t152 = -t144 * t157 + t145 * t149;
	t151 = t144 * t158 + t145 * t150;
	t143 = t154 * qJD(1) - t151 * qJD(4);
	t142 = t153 * qJD(1) + t152 * qJD(4);
	t141 = t152 * qJD(1) + t153 * qJD(4);
	t140 = t151 * qJD(1) - t154 * qJD(4);
	t1 = [0, 0, 0, -t150 * t155, 0; t141, 0, 0, t142, 0; t143, 0, 0, -t140, 0; 0, 0, 0, t149 * t155, 0; t140, 0, 0, -t143, 0; t142, 0, 0, t141, 0; 0, 0, 0, 0, 0; -t144 * t156, 0, 0, 0, 0; t145 * t156, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end