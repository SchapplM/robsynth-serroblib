% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S2RR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% Ja_rot [3x2]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:02
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S2RR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),uint8(0),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_jacobia_rot_sym_varpar: qJ has to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S2RR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_jacobia_rot_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:02:34
	% EndTime: 2019-10-09 20:02:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:02:34
	% EndTime: 2019-10-09 20:02:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN; NaN, NaN; NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:02:34
	% EndTime: 2019-10-09 20:02:34
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->2), mult. (13->3), div. (14->4), fcn. (24->3), ass. (0->5)
	t17 = sin(qJ(1));
	t19 = cos(qJ(1));
	t20 = t17 ^ 2 / t19 ^ 2;
	t11 = cos(atan2(0, t19));
	t1 = [0, 0; (0.1e1 + t20) / (0.1e1 + 0.1e1 / t11 ^ 2 * t20) / t11, 0; 0, 1;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,2);
end