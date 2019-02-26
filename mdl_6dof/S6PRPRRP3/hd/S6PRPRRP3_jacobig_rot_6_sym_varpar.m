% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRP3_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:33
% EndTime: 2019-02-26 19:51:34
% DurationCPUTime: 0.09s
% Computational Cost: add. (15->10), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->14)
t100 = sin(pkin(10));
t101 = sin(pkin(6));
t109 = t100 * t101;
t102 = cos(pkin(10));
t108 = t102 * t101;
t103 = cos(pkin(6));
t104 = sin(qJ(2));
t107 = t103 * t104;
t105 = cos(qJ(2));
t106 = t103 * t105;
t99 = pkin(11) + qJ(4);
t98 = cos(t99);
t97 = sin(t99);
t1 = [0, t109, 0, t100 * t106 + t102 * t104 (-t100 * t107 + t102 * t105) * t97 - t98 * t109, 0; 0, -t108, 0, t100 * t104 - t102 * t106 (t100 * t105 + t102 * t107) * t97 + t98 * t108, 0; 0, t103, 0, -t101 * t105, t101 * t104 * t97 - t103 * t98, 0;];
Jg_rot  = t1;
