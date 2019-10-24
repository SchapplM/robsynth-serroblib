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
% Datum: 2019-10-24 10:44
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
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
	% StartTime: 2019-10-24 10:44:26
	% EndTime: 2019-10-24 10:44:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:44:26
	% EndTime: 2019-10-24 10:44:26
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
	% StartTime: 2019-10-24 10:44:26
	% EndTime: 2019-10-24 10:44:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t43 = qJD(1) * sin(qJ(1));
	t42 = qJD(1) * cos(qJ(1));
	t39 = cos(pkin(8));
	t38 = sin(pkin(8));
	t1 = [0, 0, 0, 0, 0; t39 * t43, 0, 0, 0, 0; -t39 * t42, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t38 * t43, 0, 0, 0, 0; t38 * t42, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t42, 0, 0, 0, 0; -t43, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:44:26
	% EndTime: 2019-10-24 10:44:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (26->12), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t137 = sin(qJ(3));
	t138 = sin(qJ(1));
	t150 = t137 * t138;
	t140 = cos(qJ(1));
	t149 = t137 * t140;
	t139 = cos(qJ(3));
	t148 = t138 * t139;
	t147 = t139 * t140;
	t135 = sin(pkin(8));
	t146 = qJD(1) * t135;
	t145 = qJD(3) * t135;
	t136 = cos(pkin(8));
	t144 = t136 * t147 + t150;
	t143 = t136 * t148 - t149;
	t142 = t136 * t149 - t148;
	t141 = t136 * t150 + t147;
	t134 = t144 * qJD(1) - t141 * qJD(3);
	t133 = t142 * qJD(1) + t143 * qJD(3);
	t132 = t143 * qJD(1) + t142 * qJD(3);
	t131 = t141 * qJD(1) - t144 * qJD(3);
	t1 = [0, 0, -t139 * t145, 0, 0; t132, 0, t133, 0, 0; -t134, 0, t131, 0, 0; 0, 0, t137 * t145, 0, 0; -t131, 0, t134, 0, 0; t133, 0, t132, 0, 0; 0, 0, 0, 0, 0; t138 * t146, 0, 0, 0, 0; -t140 * t146, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:44:26
	% EndTime: 2019-10-24 10:44:26
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (119->15), mult. (132->24), div. (0->0), fcn. (132->6), ass. (0->23)
	t191 = qJD(3) + qJD(4);
	t193 = sin(pkin(8));
	t205 = t191 * t193;
	t194 = cos(pkin(8));
	t195 = sin(qJ(1));
	t204 = t194 * t195;
	t196 = cos(qJ(1));
	t203 = t194 * t196;
	t202 = qJD(1) * t193;
	t192 = qJ(3) + qJ(4);
	t190 = cos(t192);
	t201 = t190 * t205;
	t189 = sin(t192);
	t200 = t189 * t195 + t190 * t203;
	t199 = -t189 * t196 + t190 * t204;
	t198 = t189 * t203 - t190 * t195;
	t197 = t189 * t204 + t190 * t196;
	t188 = t189 * t205;
	t187 = t200 * qJD(1) - t197 * t191;
	t186 = t198 * qJD(1) + t199 * t191;
	t185 = t199 * qJD(1) + t198 * t191;
	t184 = t197 * qJD(1) - t200 * t191;
	t1 = [0, 0, -t201, -t201, 0; t185, 0, t186, t186, 0; -t187, 0, t184, t184, 0; 0, 0, t188, t188, 0; -t184, 0, t187, t187, 0; t186, 0, t185, t185, 0; 0, 0, 0, 0, 0; t195 * t202, 0, 0, 0, 0; -t196 * t202, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:44:26
	% EndTime: 2019-10-24 10:44:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (119->15), mult. (132->24), div. (0->0), fcn. (132->6), ass. (0->23)
	t199 = qJD(3) + qJD(4);
	t201 = sin(pkin(8));
	t213 = t199 * t201;
	t202 = cos(pkin(8));
	t203 = sin(qJ(1));
	t212 = t202 * t203;
	t204 = cos(qJ(1));
	t211 = t202 * t204;
	t210 = qJD(1) * t201;
	t200 = qJ(3) + qJ(4);
	t198 = cos(t200);
	t209 = t198 * t213;
	t197 = sin(t200);
	t208 = t197 * t203 + t198 * t211;
	t207 = -t197 * t204 + t198 * t212;
	t206 = t197 * t211 - t198 * t203;
	t205 = t197 * t212 + t198 * t204;
	t196 = t197 * t213;
	t195 = t208 * qJD(1) - t205 * t199;
	t194 = t206 * qJD(1) + t207 * t199;
	t193 = t207 * qJD(1) + t206 * t199;
	t192 = t205 * qJD(1) - t208 * t199;
	t1 = [0, 0, -t209, -t209, 0; t193, 0, t194, t194, 0; -t195, 0, t192, t192, 0; 0, 0, t196, t196, 0; -t192, 0, t195, t195, 0; t194, 0, t193, t193, 0; 0, 0, 0, 0, 0; t203 * t210, 0, 0, 0, 0; -t204 * t210, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end