% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR6_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:05
% EndTime: 2019-02-26 20:07:05
% DurationCPUTime: 0.04s
% Computational Cost: add. (16->14), mult. (48->33), div. (0->0), fcn. (72->10), ass. (0->18)
t126 = sin(pkin(12));
t128 = sin(pkin(6));
t140 = t126 * t128;
t127 = sin(pkin(7));
t139 = t127 * t128;
t129 = cos(pkin(12));
t138 = t129 * t128;
t131 = cos(pkin(6));
t133 = sin(qJ(2));
t137 = t131 * t133;
t135 = cos(qJ(2));
t136 = t131 * t135;
t134 = cos(qJ(3));
t132 = sin(qJ(3));
t130 = cos(pkin(7));
t125 = -t126 * t136 - t129 * t133;
t124 = -t126 * t133 + t129 * t136;
t1 = [0, t140, -t125 * t127 + t130 * t140, 0 (-t126 * t137 + t129 * t135) * t132 + (-t125 * t130 - t126 * t139) * t134, 0; 0, -t138, -t124 * t127 - t130 * t138, 0 (t126 * t135 + t129 * t137) * t132 + (-t124 * t130 + t127 * t138) * t134, 0; 0, t131, t131 * t130 - t135 * t139, 0, -t131 * t127 * t134 + (-t130 * t134 * t135 + t132 * t133) * t128, 0;];
Jg_rot  = t1;
