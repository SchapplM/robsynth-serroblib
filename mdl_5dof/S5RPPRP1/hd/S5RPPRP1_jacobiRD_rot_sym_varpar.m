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
% Datum: 2019-10-24 10:39
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
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
	% StartTime: 2019-10-24 10:39:42
	% EndTime: 2019-10-24 10:39:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:42
	% EndTime: 2019-10-24 10:39:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = qJD(1) * cos(qJ(1));
	t7 = qJD(1) * sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t7, 0, 0, 0, 0; -t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9, 0, 0, 0, 0; t7, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:42
	% EndTime: 2019-10-24 10:39:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->2), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(7);
	t13 = qJD(1) * cos(t12);
	t10 = qJD(1) * sin(t12);
	t1 = [0, 0, 0, 0, 0; t10, 0, 0, 0, 0; -t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; t13, 0, 0, 0, 0; t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:42
	% EndTime: 2019-10-24 10:39:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->5), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t49 = qJD(1) * sin(pkin(8));
	t48 = qJD(1) * cos(pkin(8));
	t45 = qJ(1) + pkin(7);
	t44 = cos(t45);
	t43 = sin(t45);
	t1 = [0, 0, 0, 0, 0; t43 * t48, 0, 0, 0, 0; -t44 * t48, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t43 * t49, 0, 0, 0, 0; t44 * t49, 0, 0, 0, 0; 0, 0, 0, 0, 0; -qJD(1) * t44, 0, 0, 0, 0; -qJD(1) * t43, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:42
	% EndTime: 2019-10-24 10:39:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (60->13), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->20)
	t142 = cos(pkin(8));
	t143 = sin(qJ(4));
	t152 = t142 * t143;
	t144 = cos(qJ(4));
	t151 = t142 * t144;
	t141 = sin(pkin(8));
	t150 = qJD(1) * t141;
	t149 = qJD(4) * t141;
	t140 = qJ(1) + pkin(7);
	t138 = sin(t140);
	t139 = cos(t140);
	t148 = t138 * t143 + t139 * t151;
	t147 = -t138 * t144 + t139 * t152;
	t146 = t138 * t151 - t139 * t143;
	t145 = t138 * t152 + t139 * t144;
	t137 = t148 * qJD(1) - t145 * qJD(4);
	t136 = t147 * qJD(1) + t146 * qJD(4);
	t135 = t146 * qJD(1) + t147 * qJD(4);
	t134 = t145 * qJD(1) - t148 * qJD(4);
	t1 = [0, 0, 0, -t144 * t149, 0; t135, 0, 0, t136, 0; -t137, 0, 0, t134, 0; 0, 0, 0, t143 * t149, 0; -t134, 0, 0, t137, 0; t136, 0, 0, t135, 0; 0, 0, 0, 0, 0; t138 * t150, 0, 0, 0, 0; -t139 * t150, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:42
	% EndTime: 2019-10-24 10:39:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (60->13), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->20)
	t152 = cos(pkin(8));
	t153 = sin(qJ(4));
	t162 = t152 * t153;
	t154 = cos(qJ(4));
	t161 = t152 * t154;
	t151 = sin(pkin(8));
	t160 = qJD(1) * t151;
	t159 = qJD(4) * t151;
	t150 = qJ(1) + pkin(7);
	t148 = sin(t150);
	t149 = cos(t150);
	t158 = t148 * t153 + t149 * t161;
	t157 = -t148 * t154 + t149 * t162;
	t156 = t148 * t161 - t149 * t153;
	t155 = t148 * t162 + t149 * t154;
	t147 = t158 * qJD(1) - t155 * qJD(4);
	t146 = t157 * qJD(1) + t156 * qJD(4);
	t145 = t156 * qJD(1) + t157 * qJD(4);
	t144 = t155 * qJD(1) - t158 * qJD(4);
	t1 = [0, 0, 0, -t154 * t159, 0; t145, 0, 0, t146, 0; -t147, 0, 0, t144, 0; 0, 0, 0, t153 * t159, 0; -t144, 0, 0, t147, 0; t146, 0, 0, t145, 0; 0, 0, 0, 0, 0; t148 * t160, 0, 0, 0, 0; -t149 * t160, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end