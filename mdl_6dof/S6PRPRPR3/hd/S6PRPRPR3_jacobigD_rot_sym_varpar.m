% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRPR3_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR3_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
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
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->4), mult. (30->12), div. (0->0), fcn. (30->8), ass. (0->10)
	t127 = sin(pkin(11));
	t130 = cos(pkin(11));
	t133 = sin(qJ(2));
	t134 = cos(qJ(2));
	t136 = qJD(2) * (t127 * t134 + t130 * t133);
	t131 = cos(pkin(10));
	t128 = sin(pkin(10));
	t126 = (t127 * t133 - t130 * t134) * qJD(2);
	t125 = cos(pkin(6)) * t136;
	t1 = [0, 0, 0, -t128 * t125 - t131 * t126, 0, 0; 0, 0, 0, t131 * t125 - t128 * t126, 0, 0; 0, 0, 0, sin(pkin(6)) * t136, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:41
	% EndTime: 2019-10-09 21:33:41
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (27->14), mult. (97->37), div. (0->0), fcn. (104->10), ass. (0->19)
	t175 = sin(pkin(11));
	t178 = cos(pkin(11));
	t182 = sin(qJ(2));
	t184 = cos(qJ(2));
	t186 = t182 * t175 - t184 * t178;
	t189 = t186 * qJD(2);
	t187 = t175 * t184 + t178 * t182;
	t173 = t187 * qJD(2);
	t177 = sin(pkin(6));
	t183 = cos(qJ(4));
	t188 = t177 * t183;
	t181 = sin(qJ(4));
	t180 = cos(pkin(6));
	t179 = cos(pkin(10));
	t176 = sin(pkin(10));
	t171 = t187 * t180;
	t170 = t180 * t189;
	t169 = t180 * t173;
	t1 = [0, 0, 0, -t176 * t169 - t179 * t189, 0, (t176 * t170 - t179 * t173) * t183 + (-(-t176 * t171 - t179 * t186) * t181 + t176 * t188) * qJD(4); 0, 0, 0, t179 * t169 - t176 * t189, 0, (-t179 * t170 - t176 * t173) * t183 + (-(t179 * t171 - t176 * t186) * t181 - t179 * t188) * qJD(4); 0, 0, 0, t177 * t173, 0, t180 * qJD(4) * t183 + (-t187 * qJD(4) * t181 - t183 * t189) * t177;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end