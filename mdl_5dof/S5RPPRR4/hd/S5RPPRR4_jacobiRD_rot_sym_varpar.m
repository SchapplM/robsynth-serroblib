% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPRR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:32:28
	% EndTime: 2020-01-03 11:32:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:32:28
	% EndTime: 2020-01-03 11:32:28
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
	% StartTime: 2020-01-03 11:32:28
	% EndTime: 2020-01-03 11:32:28
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
	% StartTime: 2020-01-03 11:32:28
	% EndTime: 2020-01-03 11:32:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t87 = cos(pkin(8));
	t88 = sin(qJ(1));
	t92 = t87 * t88;
	t89 = cos(qJ(1));
	t91 = t87 * t89;
	t90 = qJD(1) * sin(pkin(8));
	t86 = cos(pkin(9));
	t84 = sin(pkin(9));
	t1 = [0, 0, 0, 0, 0; (t84 * t89 - t86 * t92) * qJD(1), 0, 0, 0, 0; (t84 * t88 + t86 * t91) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; (t84 * t92 + t86 * t89) * qJD(1), 0, 0, 0, 0; (-t84 * t91 + t86 * t88) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; -t88 * t90, 0, 0, 0, 0; t89 * t90, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:32:28
	% EndTime: 2020-01-03 11:32:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (60->13), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->20)
	t147 = cos(pkin(8));
	t148 = sin(qJ(1));
	t157 = t147 * t148;
	t149 = cos(qJ(1));
	t156 = t147 * t149;
	t146 = sin(pkin(8));
	t155 = qJD(1) * t146;
	t154 = qJD(4) * t146;
	t145 = pkin(9) + qJ(4);
	t143 = sin(t145);
	t144 = cos(t145);
	t153 = t143 * t148 + t144 * t156;
	t152 = t143 * t149 - t144 * t157;
	t151 = -t143 * t156 + t144 * t148;
	t150 = t143 * t157 + t144 * t149;
	t142 = t153 * qJD(1) - t150 * qJD(4);
	t141 = t151 * qJD(1) + t152 * qJD(4);
	t140 = t152 * qJD(1) + t151 * qJD(4);
	t139 = t150 * qJD(1) - t153 * qJD(4);
	t1 = [0, 0, 0, -t144 * t154, 0; t140, 0, 0, t141, 0; t142, 0, 0, -t139, 0; 0, 0, 0, t143 * t154, 0; t139, 0, 0, -t142, 0; t141, 0, 0, t140, 0; 0, 0, 0, 0, 0; -t148 * t155, 0, 0, 0, 0; t149 * t155, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:32:28
	% EndTime: 2020-01-03 11:32:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (171->19), mult. (132->26), div. (0->0), fcn. (132->6), ass. (0->25)
	t197 = pkin(9) + qJ(4) + qJ(5);
	t195 = sin(t197);
	t201 = sin(qJ(1));
	t213 = t195 * t201;
	t196 = cos(t197);
	t202 = cos(qJ(1));
	t212 = t196 * t202;
	t198 = qJD(4) + qJD(5);
	t199 = sin(pkin(8));
	t211 = t198 * t199;
	t200 = cos(pkin(8));
	t210 = t200 * t201;
	t209 = t200 * t202;
	t208 = qJD(1) * t199;
	t207 = t196 * t211;
	t206 = t198 * t213;
	t205 = t198 * t212;
	t204 = t195 * t202 - t196 * t210;
	t203 = -t195 * t209 + t196 * t201;
	t192 = t195 * t211;
	t189 = -t200 * t206 - t205 + (t196 * t209 + t213) * qJD(1);
	t188 = t203 * qJD(1) + t204 * t198;
	t187 = t204 * qJD(1) + t203 * t198;
	t186 = -t200 * t205 - t206 + (t195 * t210 + t212) * qJD(1);
	t1 = [0, 0, 0, -t207, -t207; t187, 0, 0, t188, t188; t189, 0, 0, -t186, -t186; 0, 0, 0, t192, t192; t186, 0, 0, -t189, -t189; t188, 0, 0, t187, t187; 0, 0, 0, 0, 0; -t201 * t208, 0, 0, 0, 0; t202 * t208, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end