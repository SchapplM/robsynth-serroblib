% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRPR1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->4), mult. (30->12), div. (0->0), fcn. (30->8), ass. (0->10)
	t107 = sin(pkin(11));
	t110 = cos(pkin(11));
	t113 = sin(qJ(2));
	t114 = cos(qJ(2));
	t116 = qJD(2) * (t107 * t114 + t110 * t113);
	t111 = cos(pkin(10));
	t108 = sin(pkin(10));
	t106 = (t107 * t113 - t110 * t114) * qJD(2);
	t105 = cos(pkin(6)) * t116;
	t1 = [0, 0, 0, -t108 * t105 - t111 * t106, 0, 0; 0, 0, 0, t111 * t105 - t108 * t106, 0, 0; 0, 0, 0, sin(pkin(6)) * t116, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->4), mult. (30->12), div. (0->0), fcn. (30->8), ass. (0->10)
	t116 = sin(pkin(11));
	t119 = cos(pkin(11));
	t122 = sin(qJ(2));
	t123 = cos(qJ(2));
	t125 = qJD(2) * (t116 * t123 + t119 * t122);
	t120 = cos(pkin(10));
	t117 = sin(pkin(10));
	t115 = (t116 * t122 - t119 * t123) * qJD(2);
	t114 = cos(pkin(6)) * t125;
	t1 = [0, 0, 0, -t117 * t114 - t120 * t115, 0, 0; 0, 0, 0, t120 * t114 - t117 * t115, 0, 0; 0, 0, 0, sin(pkin(6)) * t125, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (36->15), mult. (97->37), div. (0->0), fcn. (104->10), ass. (0->20)
	t185 = sin(pkin(11));
	t188 = cos(pkin(11));
	t191 = sin(qJ(2));
	t192 = cos(qJ(2));
	t194 = t191 * t185 - t192 * t188;
	t197 = t194 * qJD(2);
	t195 = t185 * t192 + t188 * t191;
	t180 = t195 * qJD(2);
	t184 = qJ(4) + pkin(12);
	t182 = sin(t184);
	t187 = sin(pkin(6));
	t196 = t187 * t182;
	t190 = cos(pkin(6));
	t189 = cos(pkin(10));
	t186 = sin(pkin(10));
	t183 = cos(t184);
	t178 = t195 * t190;
	t177 = t190 * t197;
	t176 = t190 * t180;
	t1 = [0, 0, 0, -t186 * t176 - t189 * t197, 0, (t186 * t177 - t189 * t180) * t182 + ((-t186 * t178 - t189 * t194) * t183 + t186 * t196) * qJD(4); 0, 0, 0, t189 * t176 - t186 * t197, 0, (-t189 * t177 - t186 * t180) * t182 + ((t189 * t178 - t186 * t194) * t183 - t189 * t196) * qJD(4); 0, 0, 0, t187 * t180, 0, t190 * qJD(4) * t182 + (t195 * qJD(4) * t183 - t182 * t197) * t187;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end