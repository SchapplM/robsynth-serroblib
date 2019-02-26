% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRPR6_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:12
% EndTime: 2019-02-26 19:49:12
% DurationCPUTime: 0.03s
% Computational Cost: add. (8->8), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
t89 = sin(pkin(10));
t90 = sin(pkin(6));
t100 = t89 * t90;
t91 = cos(pkin(10));
t99 = t91 * t90;
t92 = cos(pkin(6));
t94 = sin(qJ(2));
t98 = t92 * t94;
t96 = cos(qJ(2));
t97 = t92 * t96;
t95 = cos(qJ(4));
t93 = sin(qJ(4));
t1 = [0, t100, 0, -t89 * t98 + t91 * t96, 0, t93 * t100 - (t89 * t97 + t91 * t94) * t95; 0, -t99, 0, t89 * t96 + t91 * t98, 0, -t93 * t99 - (t89 * t94 - t91 * t97) * t95; 0, t92, 0, t90 * t94, 0, t90 * t96 * t95 + t92 * t93;];
Jg_rot  = t1;
