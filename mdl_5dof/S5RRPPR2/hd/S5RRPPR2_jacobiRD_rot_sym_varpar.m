% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:58:03
	% EndTime: 2020-01-03 11:58:03
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:58:03
	% EndTime: 2020-01-03 11:58:03
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
	% StartTime: 2020-01-03 11:58:03
	% EndTime: 2020-01-03 11:58:03
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t23 = qJD(1) + qJD(2);
	t24 = qJ(1) + qJ(2);
	t25 = t23 * sin(t24);
	t20 = t23 * cos(t24);
	t1 = [0, 0, 0, 0, 0; -t25, -t25, 0, 0, 0; t20, t20, 0, 0, 0; 0, 0, 0, 0, 0; -t20, -t20, 0, 0, 0; -t25, -t25, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:58:03
	% EndTime: 2020-01-03 11:58:03
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (30->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t29 = qJ(1) + qJ(2) + pkin(8);
	t30 = qJD(1) + qJD(2);
	t31 = t30 * sin(t29);
	t26 = t30 * cos(t29);
	t1 = [0, 0, 0, 0, 0; -t31, -t31, 0, 0, 0; t26, t26, 0, 0, 0; 0, 0, 0, 0, 0; -t26, -t26, 0, 0, 0; -t31, -t31, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:58:03
	% EndTime: 2020-01-03 11:58:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (40->6), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t91 = qJD(1) + qJD(2);
	t97 = t91 * sin(pkin(9));
	t96 = t91 * cos(pkin(9));
	t90 = qJ(1) + qJ(2) + pkin(8);
	t88 = sin(t90);
	t95 = t88 * t96;
	t89 = cos(t90);
	t94 = t89 * t97;
	t87 = t91 * t89;
	t86 = t91 * t88;
	t85 = t89 * t96;
	t84 = t88 * t97;
	t1 = [0, 0, 0, 0, 0; -t95, -t95, 0, 0, 0; t85, t85, 0, 0, 0; 0, 0, 0, 0, 0; t84, t84, 0, 0, 0; -t94, -t94, 0, 0, 0; 0, 0, 0, 0, 0; t87, t87, 0, 0, 0; t86, t86, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:58:03
	% EndTime: 2020-01-03 11:58:03
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (171->15), mult. (132->24), div. (0->0), fcn. (132->6), ass. (0->23)
	t194 = qJD(1) + qJD(2);
	t195 = sin(pkin(9));
	t207 = t194 * t195;
	t196 = cos(pkin(9));
	t197 = sin(qJ(5));
	t206 = t196 * t197;
	t198 = cos(qJ(5));
	t205 = t196 * t198;
	t204 = qJD(5) * t195;
	t193 = qJ(1) + qJ(2) + pkin(8);
	t191 = sin(t193);
	t203 = t191 * t207;
	t192 = cos(t193);
	t202 = t191 * t197 + t192 * t205;
	t201 = t191 * t198 - t192 * t206;
	t200 = -t191 * t205 + t192 * t197;
	t199 = t191 * t206 + t192 * t198;
	t190 = t192 * t207;
	t189 = -t199 * qJD(5) + t202 * t194;
	t188 = t200 * qJD(5) + t201 * t194;
	t187 = t201 * qJD(5) + t200 * t194;
	t186 = -t202 * qJD(5) + t199 * t194;
	t1 = [0, 0, 0, 0, -t198 * t204; t187, t187, 0, 0, t188; t189, t189, 0, 0, -t186; 0, 0, 0, 0, t197 * t204; t186, t186, 0, 0, -t189; t188, t188, 0, 0, t187; 0, 0, 0, 0, 0; -t203, -t203, 0, 0, 0; t190, t190, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end