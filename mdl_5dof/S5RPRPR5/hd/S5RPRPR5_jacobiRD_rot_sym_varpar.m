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
% Datum: 2019-10-24 10:42
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
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
	% StartTime: 2019-10-24 10:42:57
	% EndTime: 2019-10-24 10:42:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:57
	% EndTime: 2019-10-24 10:42:57
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
	% StartTime: 2019-10-24 10:42:57
	% EndTime: 2019-10-24 10:42:57
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
	% StartTime: 2019-10-24 10:42:58
	% EndTime: 2019-10-24 10:42:58
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
	% StartTime: 2019-10-24 10:42:58
	% EndTime: 2019-10-24 10:42:58
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (60->13), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->20)
	t162 = cos(pkin(8));
	t163 = sin(qJ(1));
	t172 = t162 * t163;
	t164 = cos(qJ(1));
	t171 = t162 * t164;
	t161 = sin(pkin(8));
	t170 = qJD(1) * t161;
	t169 = qJD(3) * t161;
	t160 = qJ(3) + pkin(9);
	t158 = sin(t160);
	t159 = cos(t160);
	t168 = t158 * t163 + t159 * t171;
	t167 = -t158 * t164 + t159 * t172;
	t166 = t158 * t171 - t159 * t163;
	t165 = t158 * t172 + t159 * t164;
	t157 = t168 * qJD(1) - t165 * qJD(3);
	t156 = t166 * qJD(1) + t167 * qJD(3);
	t155 = t167 * qJD(1) + t166 * qJD(3);
	t154 = t165 * qJD(1) - t168 * qJD(3);
	t1 = [0, 0, -t159 * t169, 0, 0; t155, 0, t156, 0, 0; -t157, 0, t154, 0, 0; 0, 0, t158 * t169, 0, 0; -t154, 0, t157, 0, 0; t156, 0, t155, 0, 0; 0, 0, 0, 0, 0; t163 * t170, 0, 0, 0, 0; -t164 * t170, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:58
	% EndTime: 2019-10-24 10:42:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (171->15), mult. (132->24), div. (0->0), fcn. (132->6), ass. (0->23)
	t201 = qJD(3) + qJD(5);
	t202 = sin(pkin(8));
	t214 = t201 * t202;
	t203 = cos(pkin(8));
	t204 = sin(qJ(1));
	t213 = t203 * t204;
	t205 = cos(qJ(1));
	t212 = t203 * t205;
	t211 = qJD(1) * t202;
	t200 = qJ(3) + pkin(9) + qJ(5);
	t199 = cos(t200);
	t210 = t199 * t214;
	t198 = sin(t200);
	t209 = t198 * t204 + t199 * t212;
	t208 = -t198 * t205 + t199 * t213;
	t207 = t198 * t212 - t199 * t204;
	t206 = t198 * t213 + t199 * t205;
	t197 = t198 * t214;
	t196 = t209 * qJD(1) - t206 * t201;
	t195 = t207 * qJD(1) + t208 * t201;
	t194 = t208 * qJD(1) + t207 * t201;
	t193 = t206 * qJD(1) - t209 * t201;
	t1 = [0, 0, -t210, 0, -t210; t194, 0, t195, 0, t195; -t196, 0, t193, 0, t193; 0, 0, t197, 0, t197; -t193, 0, t196, 0, t196; t195, 0, t194, 0, t194; 0, 0, 0, 0, 0; t204 * t211, 0, 0, 0, 0; -t205 * t211, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end