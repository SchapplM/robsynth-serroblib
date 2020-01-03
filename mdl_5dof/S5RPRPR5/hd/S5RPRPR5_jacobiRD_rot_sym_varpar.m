% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRPR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:44:05
	% EndTime: 2020-01-03 11:44:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:44:05
	% EndTime: 2020-01-03 11:44:05
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
	% StartTime: 2020-01-03 11:44:05
	% EndTime: 2020-01-03 11:44:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t43 = qJD(1) * sin(qJ(1));
	t42 = qJD(1) * cos(qJ(1));
	t39 = cos(pkin(8));
	t38 = sin(pkin(8));
	t1 = [0, 0, 0, 0, 0; -t39 * t43, 0, 0, 0, 0; t39 * t42, 0, 0, 0, 0; 0, 0, 0, 0, 0; t38 * t43, 0, 0, 0, 0; -t38 * t42, 0, 0, 0, 0; 0, 0, 0, 0, 0; t42, 0, 0, 0, 0; t43, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:44:05
	% EndTime: 2020-01-03 11:44:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (26->12), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t133 = sin(qJ(3));
	t134 = sin(qJ(1));
	t146 = t133 * t134;
	t136 = cos(qJ(1));
	t145 = t133 * t136;
	t135 = cos(qJ(3));
	t144 = t134 * t135;
	t143 = t135 * t136;
	t131 = sin(pkin(8));
	t142 = qJD(1) * t131;
	t141 = qJD(3) * t131;
	t132 = cos(pkin(8));
	t140 = t132 * t143 + t146;
	t139 = -t132 * t144 + t145;
	t138 = -t132 * t145 + t144;
	t137 = t132 * t146 + t143;
	t130 = t140 * qJD(1) - t137 * qJD(3);
	t129 = t138 * qJD(1) + t139 * qJD(3);
	t128 = t139 * qJD(1) + t138 * qJD(3);
	t127 = t137 * qJD(1) - t140 * qJD(3);
	t1 = [0, 0, -t135 * t141, 0, 0; t128, 0, t129, 0, 0; t130, 0, -t127, 0, 0; 0, 0, t133 * t141, 0, 0; t127, 0, -t130, 0, 0; t129, 0, t128, 0, 0; 0, 0, 0, 0, 0; -t134 * t142, 0, 0, 0, 0; t136 * t142, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:44:05
	% EndTime: 2020-01-03 11:44:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (60->13), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->20)
	t158 = cos(pkin(8));
	t159 = sin(qJ(1));
	t168 = t158 * t159;
	t160 = cos(qJ(1));
	t167 = t158 * t160;
	t157 = sin(pkin(8));
	t166 = qJD(1) * t157;
	t165 = qJD(3) * t157;
	t156 = qJ(3) + pkin(9);
	t154 = sin(t156);
	t155 = cos(t156);
	t164 = t154 * t159 + t155 * t167;
	t163 = t154 * t160 - t155 * t168;
	t162 = -t154 * t167 + t155 * t159;
	t161 = t154 * t168 + t155 * t160;
	t153 = t164 * qJD(1) - t161 * qJD(3);
	t152 = t162 * qJD(1) + t163 * qJD(3);
	t151 = t163 * qJD(1) + t162 * qJD(3);
	t150 = t161 * qJD(1) - t164 * qJD(3);
	t1 = [0, 0, -t155 * t165, 0, 0; t151, 0, t152, 0, 0; t153, 0, -t150, 0, 0; 0, 0, t154 * t165, 0, 0; t150, 0, -t153, 0, 0; t152, 0, t151, 0, 0; 0, 0, 0, 0, 0; -t159 * t166, 0, 0, 0, 0; t160 * t166, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:44:05
	% EndTime: 2020-01-03 11:44:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (171->19), mult. (132->26), div. (0->0), fcn. (132->6), ass. (0->25)
	t202 = qJ(3) + pkin(9) + qJ(5);
	t200 = sin(t202);
	t206 = sin(qJ(1));
	t218 = t200 * t206;
	t201 = cos(t202);
	t207 = cos(qJ(1));
	t217 = t201 * t207;
	t203 = qJD(3) + qJD(5);
	t204 = sin(pkin(8));
	t216 = t203 * t204;
	t205 = cos(pkin(8));
	t215 = t205 * t206;
	t214 = t205 * t207;
	t213 = qJD(1) * t204;
	t212 = t201 * t216;
	t211 = t203 * t218;
	t210 = t203 * t217;
	t209 = t200 * t207 - t201 * t215;
	t208 = -t200 * t214 + t201 * t206;
	t197 = t200 * t216;
	t194 = -t205 * t211 - t210 + (t201 * t214 + t218) * qJD(1);
	t193 = t208 * qJD(1) + t209 * t203;
	t192 = t209 * qJD(1) + t208 * t203;
	t191 = -t205 * t210 - t211 + (t200 * t215 + t217) * qJD(1);
	t1 = [0, 0, -t212, 0, -t212; t192, 0, t193, 0, t193; t194, 0, -t191, 0, -t191; 0, 0, t197, 0, t197; t191, 0, -t194, 0, -t194; t193, 0, t192, 0, t192; 0, 0, 0, 0, 0; -t206 * t213, 0, 0, 0, 0; t207 * t213, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end