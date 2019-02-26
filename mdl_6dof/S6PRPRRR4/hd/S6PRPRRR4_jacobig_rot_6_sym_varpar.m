% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRR4_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:55:23
% EndTime: 2019-02-26 19:55:23
% DurationCPUTime: 0.08s
% Computational Cost: add. (26->10), mult. (39->20), div. (0->0), fcn. (63->8), ass. (0->17)
t111 = sin(pkin(11));
t112 = sin(pkin(6));
t120 = t111 * t112;
t113 = cos(pkin(11));
t119 = t113 * t112;
t114 = cos(pkin(6));
t115 = sin(qJ(2));
t118 = t114 * t115;
t116 = cos(qJ(2));
t117 = t114 * t116;
t110 = pkin(12) + qJ(4);
t109 = cos(t110);
t108 = sin(t110);
t107 = t108 * t112 * t115 - t109 * t114;
t106 = (-t111 * t118 + t113 * t116) * t108 - t109 * t120;
t105 = (t111 * t116 + t113 * t118) * t108 + t109 * t119;
t1 = [0, t120, 0, t111 * t117 + t113 * t115, t106, t106; 0, -t119, 0, t111 * t115 - t113 * t117, t105, t105; 0, t114, 0, -t112 * t116, t107, t107;];
Jg_rot  = t1;
