% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% Jg_rot [3x7]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg_rot = S7RRRRRRR1_jacobig_rot_5_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobig_rot_5_floatb_twist_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobig_rot_5_floatb_twist_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 21:20:53
% EndTime: 2018-11-26 21:20:53
% DurationCPUTime: 0.03s
% Computational Cost: add. (11->11), mult. (24->19), div. (0->0), fcn. (42->8), ass. (0->14)
t128 = sin(qJ(2));
t129 = sin(qJ(1));
t138 = t129 * t128;
t132 = cos(qJ(2));
t137 = t129 * t132;
t127 = sin(qJ(3));
t133 = cos(qJ(1));
t136 = t133 * t127;
t135 = t133 * t128;
t131 = cos(qJ(3));
t134 = t133 * t131;
t130 = cos(qJ(4));
t126 = sin(qJ(4));
t1 = [0, t129, -t135, -t129 * t131 - t132 * t136 (-t129 * t127 + t132 * t134) * t126 - t130 * t135, 0, 0; 0, -t133, -t138, -t127 * t137 + t134 (t131 * t137 + t136) * t126 - t130 * t138, 0, 0; 1, 0, t132, -t128 * t127, t128 * t131 * t126 + t132 * t130, 0, 0;];
Jg_rot  = t1;
