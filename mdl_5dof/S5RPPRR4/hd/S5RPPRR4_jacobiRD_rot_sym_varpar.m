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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
	% StartTime: 2019-12-05 17:45:55
	% EndTime: 2019-12-05 17:45:55
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:45:55
	% EndTime: 2019-12-05 17:45:55
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
	% StartTime: 2019-12-05 17:45:55
	% EndTime: 2019-12-05 17:45:55
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
	% StartTime: 2019-12-05 17:45:56
	% EndTime: 2019-12-05 17:45:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t91 = cos(pkin(8));
	t92 = sin(qJ(1));
	t96 = t91 * t92;
	t93 = cos(qJ(1));
	t95 = t91 * t93;
	t94 = qJD(1) * sin(pkin(8));
	t90 = cos(pkin(9));
	t88 = sin(pkin(9));
	t1 = [0, 0, 0, 0, 0; (-t88 * t93 + t90 * t96) * qJD(1), 0, 0, 0, 0; (-t88 * t92 - t90 * t95) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; (-t88 * t96 - t90 * t93) * qJD(1), 0, 0, 0, 0; (t88 * t95 - t90 * t92) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; t92 * t94, 0, 0, 0, 0; -t93 * t94, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:45:56
	% EndTime: 2019-12-05 17:45:56
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (60->13), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->20)
	t151 = cos(pkin(8));
	t152 = sin(qJ(1));
	t161 = t151 * t152;
	t153 = cos(qJ(1));
	t160 = t151 * t153;
	t150 = sin(pkin(8));
	t159 = qJD(1) * t150;
	t158 = qJD(4) * t150;
	t149 = pkin(9) + qJ(4);
	t147 = sin(t149);
	t148 = cos(t149);
	t157 = t147 * t152 + t148 * t160;
	t156 = -t147 * t153 + t148 * t161;
	t155 = t147 * t160 - t148 * t152;
	t154 = t147 * t161 + t148 * t153;
	t146 = t157 * qJD(1) - t154 * qJD(4);
	t145 = t155 * qJD(1) + t156 * qJD(4);
	t144 = t156 * qJD(1) + t155 * qJD(4);
	t143 = t154 * qJD(1) - t157 * qJD(4);
	t1 = [0, 0, 0, -t148 * t158, 0; t144, 0, 0, t145, 0; -t146, 0, 0, t143, 0; 0, 0, 0, t147 * t158, 0; -t143, 0, 0, t146, 0; t145, 0, 0, t144, 0; 0, 0, 0, 0, 0; t152 * t159, 0, 0, 0, 0; -t153 * t159, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:45:56
	% EndTime: 2019-12-05 17:45:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (171->15), mult. (132->24), div. (0->0), fcn. (132->6), ass. (0->23)
	t196 = qJD(4) + qJD(5);
	t197 = sin(pkin(8));
	t209 = t196 * t197;
	t198 = cos(pkin(8));
	t199 = sin(qJ(1));
	t208 = t198 * t199;
	t200 = cos(qJ(1));
	t207 = t198 * t200;
	t206 = qJD(1) * t197;
	t195 = pkin(9) + qJ(4) + qJ(5);
	t194 = cos(t195);
	t205 = t194 * t209;
	t193 = sin(t195);
	t204 = t193 * t199 + t194 * t207;
	t203 = -t193 * t200 + t194 * t208;
	t202 = t193 * t207 - t194 * t199;
	t201 = t193 * t208 + t194 * t200;
	t192 = t193 * t209;
	t191 = t204 * qJD(1) - t201 * t196;
	t190 = t202 * qJD(1) + t203 * t196;
	t189 = t203 * qJD(1) + t202 * t196;
	t188 = t201 * qJD(1) - t204 * t196;
	t1 = [0, 0, 0, -t205, -t205; t189, 0, 0, t190, t190; -t191, 0, 0, t188, t188; 0, 0, 0, t192, t192; -t188, 0, 0, t191, t191; t190, 0, 0, t189, t189; 0, 0, 0, 0, 0; t199 * t206, 0, 0, 0, 0; -t200 * t206, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end