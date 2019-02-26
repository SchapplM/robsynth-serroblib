% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRP3_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:02:23
% EndTime: 2019-02-26 20:02:23
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
t119 = sin(pkin(10));
t120 = sin(pkin(6));
t130 = t119 * t120;
t121 = cos(pkin(10));
t129 = t121 * t120;
t122 = cos(pkin(6));
t124 = sin(qJ(2));
t128 = t122 * t124;
t126 = cos(qJ(2));
t127 = t122 * t126;
t125 = cos(qJ(3));
t123 = sin(qJ(3));
t1 = [0, t130, t119 * t127 + t121 * t124, 0 (-t119 * t128 + t121 * t126) * t123 - t125 * t130, 0; 0, -t129, t119 * t124 - t121 * t127, 0 (t119 * t126 + t121 * t128) * t123 + t125 * t129, 0; 0, t122, -t120 * t126, 0, t120 * t124 * t123 - t122 * t125, 0;];
Jg_rot  = t1;
