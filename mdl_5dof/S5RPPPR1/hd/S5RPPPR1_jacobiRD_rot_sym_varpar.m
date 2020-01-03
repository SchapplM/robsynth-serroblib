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
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:21:13
	% EndTime: 2020-01-03 11:21:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:21:13
	% EndTime: 2020-01-03 11:21:13
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
	% StartTime: 2020-01-03 11:21:13
	% EndTime: 2020-01-03 11:21:13
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
	% StartTime: 2020-01-03 11:21:13
	% EndTime: 2020-01-03 11:21:13
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
	% StartTime: 2020-01-03 11:21:13
	% EndTime: 2020-01-03 11:21:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (15->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->10)
	t89 = sin(pkin(9));
	t92 = cos(pkin(8));
	t95 = t89 * t92;
	t91 = cos(pkin(9));
	t94 = t91 * t92;
	t93 = qJD(1) * sin(pkin(8));
	t88 = qJ(1) + pkin(7);
	t87 = cos(t88);
	t86 = sin(t88);
	t1 = [0, 0, 0, 0, 0; (-t86 * t94 + t87 * t89) * qJD(1), 0, 0, 0, 0; (t86 * t89 + t87 * t94) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; (t86 * t95 + t87 * t91) * qJD(1), 0, 0, 0, 0; (t86 * t91 - t87 * t95) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; -t86 * t93, 0, 0, 0, 0; t87 * t93, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:21:13
	% EndTime: 2020-01-03 11:21:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (94->14), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->21)
	t152 = qJ(1) + pkin(7);
	t148 = sin(t152);
	t154 = cos(pkin(8));
	t162 = t148 * t154;
	t150 = cos(t152);
	t161 = t150 * t154;
	t153 = sin(pkin(8));
	t160 = qJD(1) * t153;
	t159 = qJD(5) * t153;
	t151 = pkin(9) + qJ(5);
	t147 = sin(t151);
	t149 = cos(t151);
	t158 = t147 * t148 + t149 * t161;
	t157 = t147 * t150 - t149 * t162;
	t156 = -t147 * t161 + t148 * t149;
	t155 = t147 * t162 + t149 * t150;
	t146 = t158 * qJD(1) - t155 * qJD(5);
	t145 = t156 * qJD(1) + t157 * qJD(5);
	t144 = t157 * qJD(1) + t156 * qJD(5);
	t143 = t155 * qJD(1) - t158 * qJD(5);
	t1 = [0, 0, 0, 0, -t149 * t159; t144, 0, 0, 0, t145; t146, 0, 0, 0, -t143; 0, 0, 0, 0, t147 * t159; t143, 0, 0, 0, -t146; t145, 0, 0, 0, t144; 0, 0, 0, 0, 0; -t148 * t160, 0, 0, 0, 0; t150 * t160, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end