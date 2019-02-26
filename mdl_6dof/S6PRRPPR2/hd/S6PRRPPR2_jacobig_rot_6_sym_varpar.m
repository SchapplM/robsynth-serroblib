% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPPR2_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:46
% EndTime: 2019-02-26 19:58:46
% DurationCPUTime: 0.03s
% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
t96 = sin(pkin(10));
t97 = sin(pkin(6));
t107 = t96 * t97;
t98 = cos(pkin(10));
t106 = t98 * t97;
t100 = sin(qJ(2));
t105 = t96 * t100;
t101 = cos(qJ(2));
t104 = t96 * t101;
t103 = t98 * t100;
t102 = t98 * t101;
t99 = cos(pkin(6));
t95 = qJ(3) + pkin(11);
t94 = cos(t95);
t93 = sin(t95);
t1 = [0, t107, t99 * t104 + t103, 0, 0 (-t99 * t105 + t102) * t94 + t93 * t107; 0, -t106, -t99 * t102 + t105, 0, 0 (t99 * t103 + t104) * t94 - t93 * t106; 0, t99, -t97 * t101, 0, 0, t97 * t100 * t94 + t99 * t93;];
Jg_rot  = t1;
