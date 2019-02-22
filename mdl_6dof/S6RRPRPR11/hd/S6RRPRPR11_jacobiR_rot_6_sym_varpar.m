% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:28
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR11_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:28:35
% EndTime: 2019-02-22 11:28:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (85->13), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
t88 = sin(qJ(2));
t89 = sin(qJ(1));
t95 = t89 * t88;
t87 = qJ(4) + pkin(10) + qJ(6);
t85 = sin(t87);
t90 = cos(qJ(2));
t84 = t90 * t85;
t86 = cos(t87);
t94 = t90 * t86;
t91 = cos(qJ(1));
t93 = t91 * t88;
t92 = t91 * t90;
t83 = -t85 * t95 + t91 * t86;
t82 = t91 * t85 + t86 * t95;
t81 = t85 * t93 + t89 * t86;
t80 = -t89 * t85 + t86 * t93;
t1 = [t83, t85 * t92, 0, t80, 0, t80; t81, t89 * t84, 0, t82, 0, t82; 0, t88 * t85, 0, -t94, 0, -t94; -t82, t86 * t92, 0, -t81, 0, -t81; t80, t89 * t94, 0, t83, 0, t83; 0, t88 * t86, 0, t84, 0, t84; -t89 * t90, -t93, 0, 0, 0, 0; t92, -t95, 0, 0, 0, 0; 0, t90, 0, 0, 0, 0;];
JR_rot  = t1;
