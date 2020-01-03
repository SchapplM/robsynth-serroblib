% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S4PRPR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRPR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_jacobiRD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:26:01
	% EndTime: 2019-12-31 16:26:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:26:01
	% EndTime: 2019-12-31 16:26:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:26:01
	% EndTime: 2019-12-31 16:26:01
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t20 = qJD(2) * sin(qJ(2));
	t19 = qJD(2) * cos(qJ(2));
	t16 = cos(pkin(6));
	t15 = sin(pkin(6));
	t1 = [0, -t16 * t19, 0, 0; 0, -t15 * t19, 0, 0; 0, -t20, 0, 0; 0, t16 * t20, 0, 0; 0, t15 * t20, 0, 0; 0, -t19, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:26:01
	% EndTime: 2019-12-31 16:26:01
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (2->2), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t103 = qJD(2) * sin(qJ(2));
	t102 = qJD(2) * cos(qJ(2));
	t99 = cos(pkin(6));
	t98 = sin(pkin(6));
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, t99 * t102, 0, 0; 0, t98 * t102, 0, 0; 0, t103, 0, 0; 0, -t99 * t103, 0, 0; 0, -t98 * t103, 0, 0; 0, t102, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:26:02
	% EndTime: 2019-12-31 16:26:02
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (19->17), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->17)
	t148 = sin(qJ(4));
	t149 = sin(qJ(2));
	t161 = t148 * t149;
	t150 = cos(qJ(4));
	t160 = t149 * t150;
	t159 = qJD(2) * t149;
	t151 = cos(qJ(2));
	t158 = qJD(2) * t151;
	t157 = qJD(4) * t149;
	t156 = qJD(4) * t151;
	t146 = sin(pkin(6));
	t155 = t146 * t158;
	t147 = cos(pkin(6));
	t154 = t147 * t158;
	t153 = t148 * t156 + t150 * t159;
	t152 = -t148 * t159 + t150 * t156;
	t1 = [0, t152 * t147, 0, t150 * t154 + (-t146 * t150 - t147 * t161) * qJD(4); 0, t152 * t146, 0, t150 * t155 + (-t146 * t161 + t147 * t150) * qJD(4); 0, t148 * t158 + t150 * t157, 0, t153; 0, -t153 * t147, 0, -t148 * t154 + (t146 * t148 - t147 * t160) * qJD(4); 0, -t153 * t146, 0, -t148 * t155 + (-t146 * t160 - t147 * t148) * qJD(4); 0, -t148 * t157 + t150 * t158, 0, t152; 0, -t154, 0, 0; 0, -t155, 0, 0; 0, -t159, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,4);
end