% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRRP4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
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
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.03s
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
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
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
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (119->19), mult. (132->26), div. (0->0), fcn. (132->6), ass. (0->25)
	t194 = qJ(3) + qJ(4);
	t191 = sin(t194);
	t197 = sin(qJ(1));
	t209 = t191 * t197;
	t192 = cos(t194);
	t198 = cos(qJ(1));
	t208 = t192 * t198;
	t193 = qJD(3) + qJD(4);
	t195 = sin(pkin(8));
	t207 = t193 * t195;
	t196 = cos(pkin(8));
	t206 = t196 * t197;
	t205 = t196 * t198;
	t204 = qJD(1) * t195;
	t203 = t192 * t207;
	t202 = t193 * t209;
	t201 = t193 * t208;
	t200 = t191 * t198 - t192 * t206;
	t199 = -t191 * t205 + t192 * t197;
	t188 = t191 * t207;
	t185 = -t196 * t202 - t201 + (t192 * t205 + t209) * qJD(1);
	t184 = t199 * qJD(1) + t200 * t193;
	t183 = t200 * qJD(1) + t199 * t193;
	t182 = -t196 * t201 - t202 + (t191 * t206 + t208) * qJD(1);
	t1 = [0, 0, -t203, -t203, 0; t183, 0, t184, t184, 0; t185, 0, -t182, -t182, 0; 0, 0, t188, t188, 0; t182, 0, -t185, -t185, 0; t184, 0, t183, t183, 0; 0, 0, 0, 0, 0; -t197 * t204, 0, 0, 0, 0; t198 * t204, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (119->19), mult. (132->26), div. (0->0), fcn. (132->6), ass. (0->25)
	t202 = qJ(3) + qJ(4);
	t199 = sin(t202);
	t205 = sin(qJ(1));
	t217 = t199 * t205;
	t200 = cos(t202);
	t206 = cos(qJ(1));
	t216 = t200 * t206;
	t201 = qJD(3) + qJD(4);
	t203 = sin(pkin(8));
	t215 = t201 * t203;
	t204 = cos(pkin(8));
	t214 = t204 * t205;
	t213 = t204 * t206;
	t212 = qJD(1) * t203;
	t211 = t200 * t215;
	t210 = t201 * t217;
	t209 = t201 * t216;
	t208 = t199 * t206 - t200 * t214;
	t207 = -t199 * t213 + t200 * t205;
	t196 = t199 * t215;
	t193 = -t204 * t210 - t209 + (t200 * t213 + t217) * qJD(1);
	t192 = t207 * qJD(1) + t208 * t201;
	t191 = t208 * qJD(1) + t207 * t201;
	t190 = -t204 * t209 - t210 + (t199 * t214 + t216) * qJD(1);
	t1 = [0, 0, -t211, -t211, 0; t191, 0, t192, t192, 0; t193, 0, -t190, -t190, 0; 0, 0, t196, t196, 0; t190, 0, -t193, -t193, 0; t192, 0, t191, t191, 0; 0, 0, 0, 0, 0; -t205 * t212, 0, 0, 0, 0; t206 * t212, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end