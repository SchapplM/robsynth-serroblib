% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRR6_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:55
% EndTime: 2019-02-26 19:56:55
% DurationCPUTime: 0.03s
% Computational Cost: add. (13->8), mult. (39->20), div. (0->0), fcn. (63->8), ass. (0->16)
t102 = sin(pkin(11));
t103 = sin(pkin(6));
t113 = t102 * t103;
t104 = cos(pkin(11));
t112 = t104 * t103;
t105 = cos(pkin(6));
t107 = sin(qJ(2));
t111 = t105 * t107;
t109 = cos(qJ(2));
t110 = t105 * t109;
t108 = cos(qJ(4));
t106 = sin(qJ(4));
t101 = t103 * t109 * t108 + t105 * t106;
t100 = -t106 * t112 - (t102 * t107 - t104 * t110) * t108;
t99 = t106 * t113 - (t102 * t110 + t104 * t107) * t108;
t1 = [0, t113, 0, -t102 * t111 + t104 * t109, t99, t99; 0, -t112, 0, t102 * t109 + t104 * t111, t100, t100; 0, t105, 0, t103 * t107, t101, t101;];
Jg_rot  = t1;
