% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRPR6_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:37
% EndTime: 2019-02-26 21:40:37
% DurationCPUTime: 0.04s
% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
t126 = sin(pkin(6));
t131 = sin(qJ(1));
t138 = t131 * t126;
t134 = cos(qJ(1));
t137 = t134 * t126;
t125 = sin(pkin(11));
t127 = cos(pkin(11));
t130 = sin(qJ(2));
t133 = cos(qJ(2));
t136 = t133 * t125 + t130 * t127;
t135 = t130 * t125 - t133 * t127;
t132 = cos(qJ(4));
t129 = sin(qJ(4));
t128 = cos(pkin(6));
t122 = t136 * t128;
t121 = t135 * t128;
t1 = [0, t138, 0, -t131 * t121 + t134 * t136, 0 (-t131 * t122 - t134 * t135) * t132 + t129 * t138; 0, -t137, 0, t134 * t121 + t131 * t136, 0 (t134 * t122 - t131 * t135) * t132 - t129 * t137; 1, t128, 0, t135 * t126, 0, t136 * t132 * t126 + t128 * t129;];
Jg_rot  = t1;
